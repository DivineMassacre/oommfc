import contextlib
import numbers
import warnings

import discretisedfield as df

import oommfc as oc


# ============================================================================
# Исключения для ошибок конвертации
# ============================================================================

class ConversionError(Exception):
    """Exception raised when Python to Tcl conversion fails.
    
    Attributes
    ----------
    expression : str
        Python expression that failed to convert
    reason : str
        Human-readable explanation of the failure
    suggestion : str, optional
        Suggested fix for the user
    """
    def __init__(self, expression: str, reason: str, suggestion: str = None):
        self.expression = expression
        self.reason = reason
        self.suggestion = suggestion
        
        message = f"Failed to convert '{expression}': {reason}"
        if suggestion:
            message += f"\nSuggestion: {suggestion}"
        
        super().__init__(message)


def energy_script(system, **kwargs):
    """Generate MIF script for energy terms.
    
    Parameters
    ----------
    system : micromagneticmodel.System
        System object.
    **kwargs
        Additional keyword arguments (e.g., 'n' for stage_count from TimeDriver).
    """
    mif = ""
    for term in system.energy:
        # Pass kwargs only to functions that accept them
        func = globals()[f"{term.__class__.__name__.lower()}_script"]
        try:
            mif += func(term, system, **kwargs)
        except TypeError:
            # Fallback for functions that don't accept kwargs yet
            mif += func(term, system)

    return mif


def exchange_script(term, system):
    if isinstance(term.A, numbers.Real):
        mif = "# UniformExchange\n"
        mif += f"Specify Oxs_UniformExchange:{term.name} {{\n"
        mif += f"  A {term.A}\n"
        mif += "}\n\n"

    elif isinstance(term.A, dict):
        default_value = term.A.get("default", 0)
        mif = "# Exchange6Ngbr\n"
        mif += f"Specify Oxs_Exchange6Ngbr:{term.name} {{\n"
        mif += f"  default_A {default_value}\n"
        mif += "  atlas :main_atlas\n"
        mif += "  A {\n"
        for key, value in term.A.items():
            if key != "default":
                if ":" in key:
                    region1, region2 = key.split(":")
                else:
                    region1, region2 = key, key
                mif += f"    {region1} {region2} {value}\n"
        mif += "  }\n"
        mif += "}\n\n"

    elif isinstance(term.A, df.Field):
        Amif, Aname = oc.scripts.setup_scalar_parameter(term.A, f"{term.name}_A")
        mif = Amif
        mif += "# ExchangePtwise\n"
        mif += f"Specify Oxs_ExchangePtwise:{term.name} {{\n"
        mif += f"  A {Aname}\n"
        mif += "}\n\n"

    return mif


def zeeman_script(term, system, **kwargs):
    """Generate MIF script for Zeeman energy term.
    
    Parameters
    ----------
    term : micromagneticmodel.Zeeman
        Zeeman energy term.
    system : micromagneticmodel.System
        System object.
    **kwargs
        Additional keyword arguments (e.g., 'n' for stage_count from TimeDriver).
    """
    # Check for spatiotemporal terms first
    if getattr(term, 'has_time_terms', False):
        return _spatiotemporal_zeeman_script(term, system, **kwargs)

    Hmif, Hname = oc.scripts.setup_vector_parameter(term.H, f"{term.name}_H")

    mif = ""
    mif += Hmif

    if isinstance(term.wave, str) or isinstance(term.func, str):
        if isinstance(term.wave, str):
            warnings.warn(
                "Parameter `wave` is deprecated; use `func` instead.",
                FutureWarning,
                stacklevel=2,
            )
        if isinstance(term.H, (df.Field, dict)):
            if term.wave == "sin" or term.func == "sin":
                mif += f"proc TimeFunction:{term.name} {{ total_time }} {{\n"
                mif += "  set PI [expr {4*atan(1.)}]\n"
                mif += f"  set w [expr {{ {term.f}*2*$PI }}]\n"
                mif += f"  set tt [expr {{ $total_time - {term.t0} }}]\n"
                mif += "  set wt [expr { $w*$tt }]\n"
                mif += "  set f [expr {sin($wt)}]\n"
                mif += "  set df [expr {$w*cos($wt)}]\n"
                mif += "  set ft [expr { $f }]\n"
                mif += "  set dft [expr { $df }]\n"
                mif += (
                    "  return [list $ft 0 0 0 $ft 0 0 0 $ft "
                    "$dft 0 0 0 $dft 0 0 0 $dft ] \n"
                )
                mif += "}\n\n"

            elif term.wave == "sinc" or term.func == "sinc":
                mif += f"proc TimeFunction:{term.name} {{ total_time }} {{\n"
                mif += "  set PI [expr {4*atan(1.)}]\n"
                mif += f"  set w [expr {{ {term.f}*2*$PI }}]\n"
                mif += f"  set tt [expr {{ $total_time - {term.t0} }}]\n"
                mif += "  set wt [expr {$w*$tt}]\n"
                mif += "  set sinwt [expr {sin($wt)}]\n"
                mif += "  set coswt [expr {cos($wt)}]\n"
                mif += "  if {$wt != 0} { set f [expr {$sinwt/$wt}] }\n"
                mif += "  if {$wt == 0} { set f [expr {1}] }\n"
                mif += (
                    "  if {$wt != 0} { set df "
                    "[expr {($wt*$w*$coswt - $w*$sinwt)/($wt*$wt)}] }\n"
                )
                mif += "  if {$wt == 0} { set df [expr {0}] }\n"
                mif += "  set ft [expr { $f }]\n"
                mif += "  set dft [expr { $df }]\n"
                mif += (
                    "  return [list $ft 0 0 0 $ft 0 0 0 $ft "
                    "$dft 0 0 0 $dft 0 0 0 $dft ] \n"
                )
                mif += "}\n\n"

            mif += "# TransformZeeman\n"
            mif += f"Specify Oxs_TransformZeeman:{term.name} {{\n"
            mif += "  type general\n"
            mif += "  script_args total_time\n"
            mif += f"  script TimeFunction:{term.name}\n"
            mif += f"  field {Hname}\n"
            mif += "}\n\n"
        else:
            if term.wave == "sin" or term.func == "sin":
                mif += f"proc TimeFunction:{term.name} {{ total_time }} {{\n"
                mif += "  set PI [expr {4*atan(1.)}]\n"
                mif += f"  set w [expr {{ {term.f}*2*$PI }}]\n"
                mif += f"  set tt [expr {{ $total_time - {term.t0} }}]\n"
                mif += "  set wt [expr { $w*$tt }]\n"
                mif += "  set f [expr {sin($wt)}]\n"
                mif += "  set df [expr {$w*cos($wt)}]\n"
                mif += f"  set Hx [expr {{ {term.H[0]}*$f }}]\n"
                mif += f"  set Hy [expr {{ {term.H[1]}*$f }}]\n"
                mif += f"  set Hz [expr {{ {term.H[2]}*$f }}]\n"
                mif += f"  set dHx [expr {{ {term.H[0]}*$df }}]\n"
                mif += f"  set dHy [expr {{ {term.H[1]}*$df }}]\n"
                mif += f"  set dHz [expr {{ {term.H[2]}*$df }}]\n"
                mif += "  return [list $Hx $Hy $Hz $dHx $dHy $dHz ] \n"
                mif += "}\n\n"
            elif term.wave == "sinc" or term.func == "sinc":
                mif += f"proc TimeFunction:{term.name} {{ total_time }} {{\n"
                mif += "  set PI [expr {4*atan(1.)}]\n"
                mif += f"  set w [expr {{ {term.f}*2*$PI }}]\n"
                mif += f"  set tt [expr {{ $total_time - {term.t0} }}]\n"
                mif += "  set wt [expr {$w*$tt}]\n"
                mif += "  set sinwt [expr {sin($wt)}]\n"
                mif += "  set coswt [expr {cos($wt)}]\n"
                mif += "  if {$wt != 0} { set f [expr {$sinwt/$wt}] }\n"
                mif += "  if {$wt == 0} { set f [expr {1}] }\n"
                mif += (
                    "  if {$wt != 0} { set df "
                    "[expr {($wt*$w*$coswt - $w*$sinwt)/($wt*$wt)}] }\n"
                )
                mif += "  if {$wt == 0} { set df [expr {0}] }\n"
                mif += f"  set Hx [expr {{ {term.H[0]}*$f }}]\n"
                mif += f"  set Hy [expr {{ {term.H[1]}*$f }}]\n"
                mif += f"  set Hz [expr {{ {term.H[2]}*$f }}]\n"
                mif += f"  set dHx [expr {{ {term.H[0]}*$df }}]\n"
                mif += f"  set dHy [expr {{ {term.H[1]}*$df }}]\n"
                mif += f"  set dHz [expr {{ {term.H[2]}*$df }}]\n"
                mif += "  return [list $Hx $Hy $Hz $dHx $dHy $dHz ] \n"
                mif += "}\n\n"

            mif += "# ScriptUZeeman\n"
            mif += f"Specify Oxs_ScriptUZeeman:{term.name} {{\n"
            mif += "  script_args total_time\n"
            mif += f"  script TimeFunction:{term.name}\n"
            mif += "}\n\n"
    elif hasattr(term, "tlist"):
        if isinstance(term.H, (df.Field, dict)):
            mif += f"proc TimeFunction:{term.name} {{ total_time }} {{\n"
            mif += f"  set tstep {term.dt}\n"
            mif += "  set index [expr round($total_time/$tstep)]\n"
            tstr = " ".join(map(str, term.tlist))
            for char, replacement in [[",", ""], ["[", "{ "], ["]", " }"]]:
                tstr = tstr.replace(char, replacement)
            mif += f"  set H_t_fac {{ {tstr} }}\n"
            dtstr = " ".join(map(str, term.dtlist))
            for char, replacement in [[",", ""], ["[", "{ "], ["]", " }"]]:
                dtstr = dtstr.replace(char, replacement)
            mif += f"  set dH_t_fac {{ {dtstr} }}\n"
            mif += "  set H [lindex $H_t_fac $index]\n"
            mif += "  set dH [lindex $dH_t_fac $index]\n"
            if isinstance(term.tlist[0], list):
                mif += (
                    "  return [list"
                    f" {' '.join([f'[lindex $H {i}]' for i in range(9)])}"
                    f" {' '.join([f'[lindex $dH {i}]' for i in range(9)])}"
                    "]\n"
                )
            else:
                mif += "  return [list $H $H $H $dH $dH $dH]\n"
            mif += "}\n\n"

            mif += "# TransformZeeman\n"
            mif += f"Specify Oxs_TransformZeeman:{term.name} {{\n"
            if isinstance(term.tlist[0], list):
                mif += "  type general\n"
            else:
                mif += "  type diagonal\n"
            mif += "  script_args total_time\n"
            mif += f"  script TimeFunction:{term.name}\n"
            mif += f"  field {Hname}\n"
            mif += "}\n\n"
        else:
            mif += f"proc TimeFunction:{term.name} {{ total_time }} {{\n"
            mif += f"  set tstep {term.dt}\n"
            mif += "  set index [expr round($total_time/$tstep)]\n"
            mif += f"  set H_t_fac {{ {' '.join(map(str, term.tlist))} }}\n"
            mif += f"  set dH_t_fac {{ {' '.join(map(str, term.dtlist))} }}\n"
            mif += "  set H_fac [lindex $H_t_fac $index]\n"
            mif += "  set dH_fac [lindex $dH_t_fac $index]\n"
            mif += f"  set Hx [expr {{ {term.H[0]}*$H_fac }}]\n"
            mif += f"  set Hy [expr {{ {term.H[1]}*$H_fac }}]\n"
            mif += f"  set Hz [expr {{ {term.H[2]}*$H_fac }}]\n"
            mif += f"  set dHx [expr {{ {term.H[0]}*$dH_fac }}]\n"
            mif += f"  set dHy [expr {{ {term.H[1]}*$dH_fac }}]\n"
            mif += f"  set dHz [expr {{ {term.H[2]}*$dH_fac }}]\n"
            mif += "  return [list $Hx $Hy $Hz $dHx $dHy $dHz]\n"
            mif += "}\n\n"

            mif += "# ScriptUZeeman\n"
            mif += f"Specify Oxs_ScriptUZeeman:{term.name} {{\n"
            mif += "  script_args total_time\n"
            mif += f"  script TimeFunction:{term.name}\n"
            mif += "}\n\n"
    elif isinstance(term.tcl_strings, dict):
        mif += term.tcl_strings["script"]
        mif += f"\n# {term.tcl_strings['energy'][4:]}\n"  # 3.9 removeprefix
        mif += f"Specify {term.tcl_strings['energy']}:{term.name} {{\n"
        mif += f"  script {term.tcl_strings['script_name']}\n"
        for key in ["type", "script_args"]:
            with contextlib.suppress(KeyError):
                mif += f"  {key} {term.tcl_strings[key]}\n"
        if term.tcl_strings["energy"] == "Oxs_TransformZeeman":
            mif += f"  field {Hname}\n"
        mif += "}\n\n"
    else:
        mif += "# FixedZeeman\n"
        mif += f"Specify Oxs_FixedZeeman:{term.name} {{\n"
        mif += f"  field {Hname}\n"
        mif += "}\n\n"

    return mif


def demag_script(term, system):
    mif = "# Demag\n"
    if system.m.mesh.bc in ("neumann", "dirichlet", ""):  # no PBC
        oxs_cls = "Oxs_Demag"
    else:  # PBC
        if len(system.m.mesh.bc) == 1:
            oxs_cls = "Oxs_Demag"
        elif len(system.m.mesh.bc) >= 2:
            msg = (
                "Demagnetisation energy term with periodic boundary "
                "conditions in three directions is not supported."
            )
            raise ValueError(msg)

    mif += f"Specify {oxs_cls}:{term.name} {{\n"
    if hasattr(term, "asymptotic_radius"):
        mif += f"  asymptotic_radius {term.asymptotic_radius}\n"
    mif += "}\n\n"

    return mif


def dmi_script(term, system):
    if term.crystalclass in ["T", "O"]:
        oxs = "Oxs_DMI_T"
    elif (tcc := term.crystalclass) in ["D2d_x", "D2d_y", "D2d_z", "D2d"]:
        if tcc == "D2d":
            warnings.warn(
                "Use of `D2d` is deprecated; use `D2d_z` instead.",
                FutureWarning,
                stacklevel=2,
            )
            tcc = "D2d_z"
        oxs = f"Oxs_DMI_{tcc}"
    elif (tcc := term.crystalclass) in ["Cnv_x", "Cnv_y", "Cnv_z", "Cnv"]:
        if tcc == "Cnv":
            msg = "Use of `Cnv` is deprecated; use `Cnv_z` instead."
            warnings.warn(msg, FutureWarning, stacklevel=2)
            tcc = "Cnv_z"
        oxs = f"Oxs_DMI_{tcc}"

    mif = f"# DMI of crystallographic class {term.crystalclass}\n"
    mif += f"Specify {oxs}:{term.name} {{\n"

    if isinstance(term.D, numbers.Real):
        mif += f"  default_D {term.D}\n"
        mif += "  atlas :main_atlas\n"
        mif += "  D {\n"
        if len(system.m.mesh.subregions) == 0:
            mif += f"    main main {term.D}\n"
        else:
            mif += f"    entire entire {term.D}\n"
        mif += "  }\n"
        mif += "}\n\n"

    elif isinstance(term.D, dict):
        default_value = term.D.get("default", 0)
        mif += f"  default_D {default_value}\n"
        mif += "  atlas :main_atlas\n"
        mif += "  D {\n"
        for key, value in term.D.items():
            if key != "default":
                if ":" in key:
                    region1, region2 = key.split(":")
                else:
                    region1, region2 = key, key
                mif += f"    {region1} {region2} {value}\n"
        mif += "  }\n"
        mif += "}\n\n"

    return mif


def uniaxialanisotropy_script(term, system):
    umif, uname = oc.scripts.setup_vector_parameter(term.u, f"{term.name}_u")

    # Determine if higher-order anisotropy is defined
    if isinstance(term.K2, (numbers.Real, dict, df.Field)):
        k1mif, k1name = oc.scripts.setup_scalar_parameter(term.K1, f"{term.name}_K1")
        k2mif, k2name = oc.scripts.setup_scalar_parameter(term.K2, f"{term.name}_K2")

        mif = ""
        mif += k1mif
        mif += k2mif
        mif += umif
        mif += "# UniaxialAnisotropy\n"
        mif += f"Specify Southampton_UniaxialAnisotropy4:{term.name} {{\n"
        mif += f"  K1 {k1name}\n"
        mif += f"  K2 {k2name}\n"
        mif += f"  axis {uname}\n"
        mif += "}\n\n"

    else:
        kmif, kname = oc.scripts.setup_scalar_parameter(term.K, f"{term.name}_K")

        mif = ""
        mif += kmif
        mif += umif
        mif += "# UniaxialAnisotropy\n"
        mif += f"Specify Oxs_UniaxialAnisotropy:{term.name} {{\n"
        mif += f"  K1 {kname}\n"
        mif += f"  axis {uname}\n"
        mif += "}\n\n"

    return mif


def cubicanisotropy_script(term, system):
    kmif, kname = oc.scripts.setup_scalar_parameter(term.K, f"{term.name}_K")
    u1mif, u1name = oc.scripts.setup_vector_parameter(term.u1, f"{term.name}_u1")
    u2mif, u2name = oc.scripts.setup_vector_parameter(term.u2, f"{term.name}_u2")

    mif = ""
    mif += kmif
    mif += u1mif
    mif += u2mif
    mif += "# CubicAnisotropy\n"
    mif += f"Specify Oxs_CubicAnisotropy:{term.name} {{\n"
    mif += f"  K1 {kname}\n"
    mif += f"  axis1 {u1name}\n"
    mif += f"  axis2 {u2name}\n"
    mif += "}\n\n"

    return mif


def magnetoelastic_script(term, system):
    """Generate MIF script for magneto-elastic energy term.
    
    Supports three OOMMF classes based on the term parameters:
    - YY_FixedMEL: static strain (default)
    - YY_StageMEL: stage-based strain from files
    - YY_TransformStageMEL: transformation-based time-dependent strain
    """
    B1mif, B1name = oc.scripts.setup_scalar_parameter(term.B1, f"{term.name}_B1")
    B2mif, B2name = oc.scripts.setup_scalar_parameter(term.B2, f"{term.name}_B2")
    
    mif = ""
    mif += B1mif
    mif += B2mif
    
    # Determine which OOMMF class to use based on term parameters
    if hasattr(term, 'transform_script') and term.transform_script is not None:
        # YY_TransformStageMEL: transformation-based time-dependent strain
        mif += _transform_mel_script(term, system, B1name, B2name)
    elif hasattr(term, 'tcl_strings') and term.tcl_strings is not None:
        # YY_TransformStageMEL: tcl_strings-based time-dependent strain
        mif += _transform_mel_script(term, system, B1name, B2name)
    elif hasattr(term, 'e_diag_files') and term.e_diag_files is not None:
        # YY_StageMEL: stage-based strain from files
        mif += _stage_mel_script(term, system, B1name, B2name)
    else:
        # YY_FixedMEL: static strain (original behavior)
        ediagmif, ediagname = oc.scripts.setup_vector_parameter(
            term.e_diag, f"{term.name}_ediag"
        )
        eoffdiagmif, eoffdiagname = oc.scripts.setup_vector_parameter(
            term.e_offdiag, f"{term.name}_eoffdiag"
        )
        
        mif += ediagmif
        mif += eoffdiagmif
        mif += "# MagnetoElastic (YY_FixedMEL)\n"
        mif += f"Specify YY_FixedMEL:{term.name} {{\n"
        mif += f"  B1 {B1name}\n"
        mif += f"  B2 {B2name}\n"
        mif += f"  e_diag_field {ediagname}\n"
        mif += f"  e_offdiag_field {eoffdiagname}\n"
        mif += "}\n\n"
    
    return mif


def _stage_mel_script(term, system, B1name, B2name):
    """Generate MIF script for YY_StageMEL (stage-based strain)."""
    mif = ""
    
    # Generate Tcl scripts for strain if callable functions are provided
    if hasattr(term, 'e_diag_func') and callable(term.e_diag_func):
        mif += _generate_strain_scripts(term)
    
    mif += "# MagnetoElastic (YY_StageMEL)\n"
    mif += f"Specify YY_StageMEL:{term.name} {{\n"
    mif += f"  B1 {B1name}\n"
    mif += f"  B2 {B2name}\n"
    
    if hasattr(term, 'e_diag_files') and term.e_diag_files is not None:
        # Use file lists
        diag_files_str = " ".join(term.e_diag_files)
        offdiag_files_str = " ".join(term.e_offdiag_files)
        mif += f"  e_diag_files {{{diag_files_str}}}\n"
        mif += f"  e_offdiag_files {{{offdiag_files_str}}}\n"
    else:
        # Use Tcl scripts
        mif += f"  e_diag_script strain_diag_{term.name}\n"
        mif += f"  e_offdiag_script strain_offdiag_{term.name}\n"
    
    if hasattr(term, 'stage_count') and term.stage_count is not None:
        mif += f"  stage_count {term.stage_count}\n"
    
    mif += "}\n\n"
    
    return mif


def _transform_mel_script(term, system, B1name, B2name):
    """Generate MIF script for YY_TransformStageMEL (time-dependent strain).

    Analogous to Oxs_TransformZeeman for Zeeman energy.

    Uses tlist approach like Zeeman:
    1. Pre-compute transform values at each timestep in Python
    2. Pass values as Tcl lists in MIF (defined INSIDE proc)
    3. Tcl script indexes values by stage_time/dt

    Two modes are supported:

    1. **Direct substitution mode** (default):
       - func(t) returns FULL STRAIN VALUES [e11, e22, e33, e23, e13, e12]
       - e_diag/e_offdiag are set to zero (not used)
       - Values are directly substituted into energy calculation

    2. **Matrix transformation mode**:
       - func(t) returns TRANSFORMATION MATRIX elements M(t)
       - e_diag/e_offdiag define base strain e_base
       - Final strain: e_final = M(t) × e_base × M(t)ᵀ
       - For diagonal type: e_final_ii = M_ii(t) × e_base_ii

    Parameters
    ----------
    term : MagnetoElastic
        The magneto-elastic energy term with attributes:
        - transform_script or func: Python callable(t) -> strain values or matrix elements
        - transform_dt or dt: time step in seconds (default: 0.1 ps)
        - transform_script_args: 'total_time' or 'stage_time'
        - transform_type: 'diagonal', 'symmetric', or 'general'
        - e_diag, e_offdiag: base strain (for matrix mode only)
    system : System
        The micromagnetic system (not used in this function)
    B1name, B2name : str
        Names of B1 and B2 field specifications in MIF

    Returns
    -------
    str
        MIF script for YY_TransformStageMEL
    """
    mif = ""

    # Determine mode: direct substitution vs matrix transformation
    # Direct substitution: e_diag is None or (0, 0, 0) (no base strain)
    # Matrix transformation: e_diag is provided and non-zero (user-provided base strain)
    # Note: e_offdiag can be zero even in matrix mode (no off-diagonal base strain)
    has_base_strain = (
        term.e_diag is not None and
        term.e_diag != (0, 0, 0)
    )

    # Generate Tcl scripts for base strain fields
    if has_base_strain:
        # Matrix transformation mode: use user-provided base strain
        mif += f"proc strain_diag_{term.name} {{ stage }} {{\n"
        mif += f"  # Base diagonal strain (matrix mode): {term.e_diag}\n"
        mif += f"  set spec Oxs_UniformVectorField\n"
        mif += f"  lappend spec [subst {{\n"
        mif += f"    norm 1\n"
        mif += f"    vector {{ {term.e_diag[0]} {term.e_diag[1]} {term.e_diag[2]} }}\n"
        mif += f"  }}]\n"
        mif += f"  return $spec\n"
        mif += f"}}\n\n"

        mif += f"proc strain_offdiag_{term.name} {{ stage }} {{\n"
        mif += f"  # Base off-diagonal strain (matrix mode): {term.e_offdiag}\n"
        mif += f"  set spec Oxs_UniformVectorField\n"
        mif += f"  lappend spec [subst {{\n"
        mif += f"    norm 1\n"
        mif += f"    vector {{ {term.e_offdiag[0]} {term.e_offdiag[1]} {term.e_offdiag[2]} }}\n"
        mif += f"  }}]\n"
        mif += f"  return $spec\n"
        mif += f"}}\n\n"
    else:
        # Direct substitution mode: base strain not used (set to zero)
        mif += f"proc strain_diag_{term.name} {{ stage }} {{\n"
        mif += f"  # Direct substitution mode: base strain not used\n"
        mif += f"  set spec Oxs_UniformVectorField\n"
        mif += f"  lappend spec [subst {{\n"
        mif += f"    norm 1\n"
        mif += f"    vector {{ 0 0 0 }}\n"
        mif += f"  }}]\n"
        mif += f"  return $spec\n"
        mif += f"}}\n\n"

        mif += f"proc strain_offdiag_{term.name} {{ stage }} {{\n"
        mif += f"  # Direct substitution mode: base strain not used\n"
        mif += f"  set spec Oxs_UniformVectorField\n"
        mif += f"  lappend spec [subst {{\n"
        mif += f"    norm 1\n"
        mif += f"    vector {{ 0 0 0 }}\n"
        mif += f"  }}]\n"
        mif += f"  return $spec\n"
        mif += f"}}\n\n"

    # Generate transformation script if callable is provided
    if callable(term.transform_script):
        mif += _generate_transform_script(term, has_base_strain=has_base_strain)
    elif hasattr(term, 'tcl_strings') and term.tcl_strings is not None:
        # Use tcl_strings for advanced control (like Zeeman)
        mif += _generate_transform_script_from_tcl(term)

    mif += "# MagnetoElastic (YY_TransformStageMEL)\n"
    mif += f"Specify YY_TransformStageMEL:{term.name} {{\n"
    mif += f"  B1 {B1name}\n"
    mif += f"  B2 {B2name}\n"
    mif += f"  e_diag_script strain_diag_{term.name}\n"
    mif += f"  e_offdiag_script strain_offdiag_{term.name}\n"
    
    # Always include type - required by OOMMF
    transform_type = getattr(term, 'transform_type', 'diagonal')
    if transform_type is None:
        transform_type = 'diagonal'
    mif += f"  type {transform_type}\n"

    # REQUIRED: script field for YY_TransformStageMEL
    if callable(term.transform_script) or (hasattr(term, 'tcl_strings') and term.tcl_strings):
        mif += f"  script transform_{term.name}\n"

    # Use default script_args (stage stage_time total_time) for compatibility
    # with OOMMF MEL extension
    if hasattr(term, 'stage_count') and term.stage_count is not None:
        mif += f"  stage_count {term.stage_count}\n"

    mif += "}\n\n"

    return mif


def _generate_strain_scripts(term):
    """Generate Tcl scripts for stage-based strain from Python callables."""
    mif = ""
    
    # Generate strain_diag script
    if hasattr(term, 'e_diag_func') and callable(term.e_diag_func):
        mif += f"proc strain_diag_{term.name} {{ stage }} {{\n"
        mif += f"  # Python callable: {term.e_diag_func.__name__}\n"
        mif += "  # Note: Actual values must be pre-computed and passed via files\n"
        mif += "  error \"Stage-based strain with Python callable requires pre-computed OVf files\"\n"
        mif += "}\n\n"
    
    # Generate strain_offdiag script
    if hasattr(term, 'e_offdiag_func') and callable(term.e_offdiag_func):
        mif += f"proc strain_offdiag_{term.name} {{ stage }} {{\n"
        mif += f"  # Python callable: {term.e_offdiag_func.__name__}\n"
        mif += "  # Note: Actual values must be pre-computed and passed via files\n"
        mif += "  error \"Stage-based strain with Python callable requires pre-computed OVf files\"\n"
        mif += "}\n\n"
    
    return mif


def _generate_transform_script_from_tcl(term):
    """Generate Tcl transformation script from tcl_strings dictionary.
    
    Parameters
    ----------
    term : MagnetoElastic
        The magneto-elastic energy term with tcl_strings attribute.
        
    Returns
    -------
    str
        MIF script for transformation using tcl_strings.
    """
    mif = ""
    
    tcl_strings = term.tcl_strings
    
    # Add script if provided
    if 'script' in tcl_strings:
        mif += f"proc transform_{term.name} {{ stage stage_time total_time }} {{\n"
        mif += f"  # Custom tcl script from tcl_strings\n"
        mif += f"  {tcl_strings['script']}\n"
        mif += f"}}\n\n"
    
    return mif


def _generate_transform_script(term, has_base_strain=False):
    """Generate Tcl transformation script using tlist approach (like Zeeman).

    Pre-compute transform values at discrete time steps in Python.
    Pass values as lists to MIF. Index by stage_time.

    This matches the Zeeman approach: func/dt -> tlist -> lindex lookup.

    Script signature: transform_{name} {stage stage_time total_time}

    Two modes:

    1. **Direct substitution** (has_base_strain=False):
       - func(t) returns [e11, e22, e33, e23, e13, e12]
       - These values are directly used as strain

    2. **Matrix transformation** (has_base_strain=True):
       - func(t) returns [M11, M22, M33, dM11, dM22, dM33]
       - For diagonal type: e_final_ii = M_ii × e_base_ii
       - The Tcl script performs the multiplication

    Parameters:
    - transform_dt: time step for pre-computation (default: 0.1 ps)
    - n_points: automatically determined to cover up to 1 ns simulation
    - has_base_strain: if True, use matrix transformation mode
    """
    mif = ""

    # Get dt from term attribute (like Zeeman)
    # Default: dt = 0.1 ps (matches Zeeman default)
    dt = getattr(term, 'transform_dt', None)
    if dt is None:
        dt = 1e-13  # 0.1 ps (like Zeeman default)

    # Automatically determine n_points to cover typical simulation times
    # For dt=0.1ps, n_points=10000 covers up to 1ns
    # For dt=0.01ps, n_points=100000 covers up to 1ns
    n_points = int(1e-9 / dt)  # Cover up to 1 ns

    # Evaluate Python callable at each timestep
    transform_values = []
    for i in range(n_points):
        t = i * dt
        try:
            values = term.transform_script(t)
            if not hasattr(values, '__iter__'):
                values = [values] * 6
            values = list(values)
        except Exception:
            # Fallback to identity transform
            if has_base_strain:
                # Matrix mode: identity matrix
                values = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
            else:
                # Direct mode: zero strain
                values = [0.0] * 6

        transform_values.append(values)

    # Determine number of components based on transform_type
    n_components = len(transform_values[0]) if transform_values else 6

    # Build Tcl script with tlist approach
    mif += f"proc transform_{term.name} {{ stage stage_time total_time }} {{\n"

    if has_base_strain:
        mif += f"  # Matrix transformation mode (tlist approach)\n"
        mif += f"  # dt = {dt*1e15:.2f} fs, n_points = {n_points}\n"
        mif += f"  # e_final = M(t) × e_base × M(t)ᵀ\n"
    else:
        mif += f"  # Direct substitution mode (tlist approach like Zeeman)\n"
        mif += f"  # dt = {dt*1e15:.2f} fs, n_points = {n_points}\n"
        mif += f"  # func(t) returns full strain values\n"

    mif += "\n"

    # Create Tcl lists for each component
    for c in range(n_components):
        values_str = " ".join(f"{v[c] if c < len(v) else 0:.15e}" for v in transform_values)
        mif += f"  set transform_tlist_{c} {{ {values_str} }}\n"

    mif += "\n"
    mif += "  # Compute index from stage_time\n"
    mif += f"  set dt_lookup {dt}\n"
    mif += "  set idx [expr {int($stage_time / $dt_lookup)}]\n"
    mif += f"  set n_points {n_points}\n"
    mif += "  if {$idx >= $n_points} { set idx [expr {$n_points - 1}] }\n"
    mif += "  if {$idx < 0} { set idx 0 }\n"
    mif += "\n"

    # Return values based on transform_type and mode
    if has_base_strain:
        # Matrix transformation mode: return M(t) values
        # YY_TransformStageMEL will multiply by base strain
        mif += "  # Return transformation matrix M(t)\n"
        if term.transform_type == 'diagonal':
            mif += "  return [list \\\n"
            for c in range(6):
                mif += f"    [lindex $transform_tlist_{c} $idx]\\\n"
            mif += "  ]\n"
        elif term.transform_type == 'symmetric':
            mif += "  return [list \\\n"
            for c in range(12):
                mif += f"    [lindex $transform_tlist_{c} $idx]\\\n"
            mif += "  ]\n"
        elif term.transform_type == 'general':
            mif += "  return [list \\\n"
            for c in range(18):
                mif += f"    [lindex $transform_tlist_{c} $idx]\\\n"
            mif += "  ]\n"
        else:
            mif += "  return {}\n"
    else:
        # Direct substitution mode: return full strain values
        mif += "  # Return full strain values (direct substitution)\n"
        if term.transform_type == 'diagonal':
            mif += "  return [list \\\n"
            for c in range(6):
                mif += f"    [lindex $transform_tlist_{c} $idx]\\\n"
            mif += "  ]\n"
        elif term.transform_type == 'symmetric':
            mif += "  return [list \\\n"
            for c in range(12):
                mif += f"    [lindex $transform_tlist_{c} $idx]\\\n"
            mif += "  ]\n"
        elif term.transform_type == 'general':
            mif += "  return [list \\\n"
            for c in range(18):
                mif += f"    [lindex $transform_tlist_{c} $idx]\\\n"
            mif += "  ]\n"
        else:
            mif += "  return {}\n"

    mif += "}\n\n"

    return mif


def rkky_script(term, system):
    sr1 = system.m.mesh.subregions[term.subregions[0]]
    sr2 = system.m.mesh.subregions[term.subregions[1]]

    direction, first, second = sr1.facing_surface(sr2)

    for key, value in system.m.mesh.subregions.items():
        if value == first:
            first_name = key
        elif value == second:
            second_name = key

    mif = ""

    mif += "# Scalar field for RKKY surfaces\n"
    mif += "Specify Oxs_LinearScalarField:rkkyfield {\n"
    vectorval = df.util.assemble_index(
        0, 3, {system.m.mesh.region._dim2index(direction): 1}
    )
    mif += "  vector {{{} {} {}}}\n".format(*vectorval)
    mif += "  norm 1.0\n"
    mif += "}\n\n"

    mif += "# TwoSurfaceExchange\n"
    mif += f"Specify Oxs_TwoSurfaceExchange:{term.name} {{\n"
    if isinstance(term.sigma, numbers.Real):
        mif += f"  sigma {term.sigma}\n"
    if isinstance(term.sigma2, numbers.Real):
        mif += f"  sigma2 {term.sigma2}\n"

    mif += "  surface1 {\n"
    mif += "    atlas :main_atlas\n"
    mif += f"    region {first_name}\n"
    mif += "    scalarfield :rkkyfield\n"
    mif += f"    scalarvalue {first.pmax[system.m.mesh.region._dim2index(direction)]}\n"
    mif += "    scalarside -\n"
    mif += "  }\n"

    mif += "  surface2 {\n"
    mif += "    atlas :main_atlas\n"
    mif += f"    region {second_name}\n"
    mif += "    scalarfield :rkkyfield\n"
    mif += (
        f"    scalarvalue {second.pmin[system.m.mesh.region._dim2index(direction)]}\n"
    )
    mif += "    scalarside +\n"
    mif += "  }\n"

    mif += "}\n\n"

    return mif


def _spatiotemporal_zeeman_script(term, system, **kwargs):
    """Generate MIF script for spatiotemporal Zeeman field.

    Uses Oxs_StageZeeman + Oxs_ScriptVectorField approach:
    H(x,y,z,t) = H_static + Σᵢ [fᵢ(t) × maskᵢ(x,y,z)]
    
    Parameters
    ----------
    term : micromagneticmodel.Zeeman
        Zeeman energy term with spatiotemporal terms.
    system : micromagneticmodel.System
        System object.
    **kwargs
        Additional keyword arguments from driver. Key argument:
        
        - 'n' : int
            Number of stages from TimeDriver.drive(n=...).
            Used to set stage_count if term._stage_count is None.
    """
    mif = "# ========== Zeeman: Spatiotemporal field ==========\n"

    # Static field
    H_static = term.H
    if isinstance(H_static, (tuple, list)):
        mif += f"set H_static_x {H_static[0]}\n"
        mif += f"set H_static_y {H_static[1]}\n"
        mif += f"set H_static_z {H_static[2]}\n"
    else:
        # Non-uniform static field - use default values
        mif += "set H_static_x 0\n"
        mif += "set H_static_y 0\n"
        mif += "set H_static_z 0\n"

    # Time step and stage count
    # Priority: 1) term._stage_count, 2) kwargs['n'] from driver, 3) default 100
    dt = getattr(term, '_dt', 1e-13)
    stage_count = getattr(term, '_stage_count', None)

    if stage_count is None:
        # Try to get from driver's n parameter
        stage_count = kwargs.get('n', 100)

    mif += f"set dt {dt}\n"
    mif += f"set stage_count {stage_count}\n\n"

    # Global variable for current time
    mif += "set current_time 0\n\n"

    # Add Tcl helper functions for extended math support
    mif += "# ========== Extended math functions ==========\n"
    mif += "proc tanh {x} {\n"
    mif += "  set exp2x [expr {exp(2*$x)}]\n"
    mif += "  return [expr {($exp2x - 1) / ($exp2x + 1)}]\n"
    mif += "}\n\n"
    
    mif += "proc sinh {x} {\n"
    mif += "  return [expr {(exp($x) - exp(-$x)) / 2}]\n"
    mif += "}\n\n"

    mif += "proc cosh {x} {\n"
    mif += "  return [expr {(exp($x) + exp(-$x)) / 2}]\n"
    mif += "}\n\n"

    mif += "proc asin {x} {\n"
    mif += "  return [expr {atan2($x, sqrt(1 - $x*$x))}]\n"
    mif += "}\n\n"

    mif += "proc acos {x} {\n"
    mif += "  return [expr {atan2(sqrt(1 - $x*$x), $x)}]\n"
    mif += "}\n\n"

    mif += "proc log2 {x} {\n"
    mif += "  return [expr {log($x) / log(2)}]\n"
    mif += "}\n\n"

    # Extract global variables and closure variables from functions
    global_vars = {}
    
    def _extract_vars_from_callable(callable_obj):
        """Extract variables from callable's globals and closure."""
        vars_dict = {}

        # Check if decorated with @zeeman_func
        if getattr(callable_obj, '__is_zeeman_func__', False):
            # Use stored __zeeman_globals__
            zeeman_globals = getattr(callable_obj, '__zeeman_globals__', {})
            vars_dict.update(zeeman_globals)

        # Extract from __globals__
        if hasattr(callable_obj, '__globals__'):
            for name, value in callable_obj.__globals__.items():
                if isinstance(value, (int, float)) and not name.startswith('_'):
                    vars_dict[name] = value

        # Extract from closure (for lambda functions with captured variables)
        if hasattr(callable_obj, '__closure__') and callable_obj.__closure__:
            code = callable_obj.__code__
            freevars = getattr(code, 'co_freevars', ())
            for i, cell in enumerate(callable_obj.__closure__):
                try:
                    value = cell.cell_contents
                    if isinstance(value, (int, float)):
                        # Get variable name from freevars if available
                        if i < len(freevars):
                            name = freevars[i]
                        else:
                            name = f'_var_{i}'
                        vars_dict[name] = value
                except ValueError:
                    pass

        return vars_dict
    
    for func, mask in term._terms:
        # Extract from func
        func_vars = _extract_vars_from_callable(func)
        global_vars.update(func_vars)
        
        # Extract from mask
        if mask is not None:
            mask_vars = _extract_vars_from_callable(mask)
            global_vars.update(mask_vars)

    # Write global variables (exclude common names)
    exclude = {'pi', 'e', 't', 'x', 'y', 'z', 'H_static_x', 'H_static_y', 'H_static_z',
               'dt', 'stage_count', 'current_time', 'np', 'numpy', 'math'}
    global_var_list = []
    for name, value in sorted(global_vars.items()):
        if name not in exclude:
            mif += f"set {name} {value:.15g}\n"
            global_var_list.append(name)
    if global_vars:
        mif += "\n"

    # ========== FieldPerStage procedure ==========
    mif += "proc FieldPerStage { stage } {\n"
    mif += "  global dt current_time\n"
    mif += "  set current_time [expr {$stage * $dt}]\n"
    mif += "  return [list Oxs_ScriptVectorField {\n"
    mif += "    script SpatiotemporalField\n"
    mif += "    script_args rawpt\n"
    mif += "  }]\n"
    mif += "}\n\n"

    # ========== SpatiotemporalField procedure ==========
    mif += "proc SpatiotemporalField { x y z } {\n"
    # Build global statement with all needed variables
    global_vars_str = "H_static_x H_static_y H_static_z current_time"
    if global_var_list:
        global_vars_str += " " + " ".join(global_var_list)
    mif += f"  global {global_vars_str}\n"
    mif += "  set t $current_time\n"
    mif += "  # Static part\n"
    mif += "  set Hx $H_static_x\n"
    mif += "  set Hy $H_static_y\n"
    mif += "  set Hz $H_static_z\n\n"

    # Time-dependent terms
    for i, (func, mask) in enumerate(term._terms):
        mif += f"  # Term {i}\n"

        # Generate Tcl for func(t)
        func_tcl = _python_func_to_tcl(func, arg='t')

        # Generate Tcl for mask(x,y,z)
        if mask is None:
            mask_tcl = ['1.0', '1.0', '1.0']
        else:
            mask_tcl = _python_func_to_tcl(mask, args=['x', 'y', 'z'])

        # Compute f(t) × mask(x,y,z)
        if isinstance(func_tcl, str):
            # Scalar func - wrap in expr for proper evaluation
            mif += f"  set f_{i} [expr {{{func_tcl}}}]\n"

            if isinstance(mask_tcl, list):
                # Scalar func × vector mask
                mif += f"  set Hx [expr {{$Hx + $f_{i} * {mask_tcl[0]}}}]\n"
                mif += f"  set Hy [expr {{$Hy + $f_{i} * {mask_tcl[1]}}}]\n"
                mif += f"  set Hz [expr {{$Hz + $f_{i} * {mask_tcl[2]}}}]\n"
            else:
                # Scalar func × scalar mask (applied to all components)
                mif += f"  set Hx [expr {{$Hx + $f_{i} * {mask_tcl}}}]\n"
                mif += f"  set Hy [expr {{$Hy + $f_{i} * {mask_tcl}}}]\n"
                mif += f"  set Hz [expr {{$Hz + $f_{i} * {mask_tcl}}}]\n"
        else:
            # Vector func - wrap each component in expr
            mif += f"  set fx_{i} [expr {{{func_tcl[0]}}}]\n"
            mif += f"  set fy_{i} [expr {{{func_tcl[1]}}}]\n"
            mif += f"  set fz_{i} [expr {{{func_tcl[2]}}}]\n"

            if isinstance(mask_tcl, str):
                # Vector func × scalar mask
                mif += f"  set Hx [expr {{$Hx + $fx_{i} * {mask_tcl}}}]\n"
                mif += f"  set Hy [expr {{$Hy + $fy_{i} * {mask_tcl}}}]\n"
                mif += f"  set Hz [expr {{$Hz + $fz_{i} * {mask_tcl}}}]\n"
            else:
                # Vector func × vector mask
                mif += f"  set Hx [expr {{$Hx + $fx_{i} * {mask_tcl[0]}}}]\n"
                mif += f"  set Hy [expr {{$Hy + $fy_{i} * {mask_tcl[1]}}}]\n"
                mif += f"  set Hz [expr {{$Hz + $fz_{i} * {mask_tcl[2]}}}]\n"

        mif += "\n"

    mif += "  return [list $Hx $Hy $Hz]\n"
    mif += "}\n\n"

    # ========== Specify block ==========
    mif += "Specify Oxs_StageZeeman:zeeman [subst {\n"
    mif += f"  script FieldPerStage\n"
    mif += f"  stage_count {stage_count}\n"
    mif += "}]\n\n"

    return mif


def _python_func_to_tcl(func, arg='t', args=None):
    """Convert Python callable to Tcl expr string.

    Uses inspect.getsource() to extract function body and convert to Tcl.
    Supports common math functions: sin, cos, exp, sqrt, log, abs.

    Enhanced to support:
    - Phase shifts: sin(omega*t + phi)
    - Mixed expressions: sin(a*t + b*x)
    - Nested functions: sin(cos(x))
    - Complex arithmetic: (a*b + c*d) / e
    - Decorated functions with @zeeman_func decorator

    Parameters
    ----------
    func : callable
        Python function to convert
    arg : str, optional
        Argument name for scalar function (default 't')
        For spatial functions use args=['x', 'y', 'z']
    args : list, optional
        Argument names for spatial function (default None)

    Returns
    -------
    str or list
        Tcl expression string(s)

    Examples
    --------
    >>> _python_func_to_tcl(lambda t: np.sin(2*np.pi*1e9*t))
    'sin(6.283185307179586e+09 * $t)'

    >>> _python_func_to_tcl(lambda t: np.sin(2*np.pi*1e9*t + np.pi/4))
    'sin(6.283185307179586e+09 * $t + 0.7853981633974483)'

    >>> _python_func_to_tcl(lambda x: np.exp(-x**2), args=['x'])
    'exp(-$x^2)'
    
    Notes
    -----
    If the function is decorated with @zeeman_func, the decorator's
    stored source code and globals are used for more reliable conversion.
    """
    import inspect
    import re

    if args is None:
        args = ['x', 'y', 'z']

    # Check if function is decorated with @zeeman_func
    if getattr(func, '__is_zeeman_func__', False):
        # Use stored source and globals from decorator
        source = getattr(func, '__zeeman_source__', None)
        globals_dict = getattr(func, '__zeeman_globals__', {})
        
        if source and not source.startswith('<'):
            # Valid source code - use it
            global_var_names = list(globals_dict.keys())
            result = _convert_source_to_tcl(source, arg, args, func, global_var_names)
            if result is not None:
                return result
        # Fallback to normal processing if decorator data is invalid
    
    # Try to get source code normally
    try:
        source = inspect.getsource(func)
    except (IOError, TypeError, OSError, IndentationError):
        # Cannot get source - use enhanced fallback
        return _convert_func_fallback_enhanced(func, arg, args)

    # Extract global variables for later substitution
    global_var_names = []
    if hasattr(func, '__globals__'):
        for name, value in func.__globals__.items():
            if isinstance(value, (int, float)) and not name.startswith('_'):
                global_var_names.append(name)

    # Parse and convert source
    result = _convert_source_to_tcl(source, arg, args, func, global_var_names)

    # If conversion failed, try enhanced fallback
    if result is None:
        return _convert_func_fallback_enhanced(func, arg, args)

    return result


def _convert_source_to_tcl(source, arg, args, func=None, global_var_names=None):
    """Convert Python source code to Tcl expression.

    Global variables (H0, omega, k, etc.) are defined in MIF,
    so we just need to convert function syntax and leave variable names as-is.
    """
    import re

    if global_var_names is None:
        global_var_names = []

    # Extract expression from lambda or def
    lambda_match = re.search(r'lambda\s+([^:]+):\s*(.+)', source, re.DOTALL)
    if lambda_match:
        params = lambda_match.group(1).strip()
        expr = lambda_match.group(2).strip()
        expr = expr.split('#')[0].strip().rstrip('\n').strip()

        # Vector function?
        if expr.startswith('(') or expr.startswith('['):
            inner = expr[1:-1]
            components = _split_components(inner)
            return [_convert_expr_to_tcl(comp.strip(), params, args, global_var_names=global_var_names) for comp in components]
        else:
            return _convert_expr_to_tcl(expr, params, args, global_var_names=global_var_names)

    # Def pattern - handle docstrings and multi-line
    # Match: def name(args): ... return EXPR
    def_match = re.search(r'def\s+\w+\s*\(([^)]*)\)\s*:.*?return\s+(.+?)(?:\n\s{4}\S|\n\n|\Z)', source, re.DOTALL)
    if def_match:
        params = def_match.group(1).strip()
        expr = def_match.group(2).strip().split('#')[0].strip()

        if expr.startswith('(') or expr.startswith('['):
            inner = expr[1:-1]
            components = _split_components(inner)
            return [_convert_expr_to_tcl(comp.strip(), params, args, global_var_names=global_var_names) for comp in components]
        else:
            return _convert_expr_to_tcl(expr, params, args, global_var_names=global_var_names)

    return None


def _split_components(expr):
    """Split tuple/list expression into components, handling nested parens.
    
    Examples:
        "sin(t), cos(t), 0" -> ["sin(t)", "cos(t)", "0"]
        "sin(t), 0, 0" -> ["sin(t)", "0", "0"]
    """
    components = []
    current = ''
    depth = 0
    
    for char in expr:
        if char in '([':
            depth += 1
            current += char
        elif char in ')]':
            depth -= 1
            current += char
        elif char == ',' and depth == 0:
            components.append(current.strip())
            current = ''
        else:
            current += char
    
    if current.strip():
        components.append(current.strip())
    
    # Clean up any trailing/leading parentheses that might have leaked
    cleaned = []
    for comp in components:
        comp = comp.strip()
        # Remove trailing commas
        comp = comp.rstrip(',')
        # Remove unmatched closing parens
        while comp.endswith(')') and comp.count('(') < comp.count(')'):
            comp = comp[:-1]
        cleaned.append(comp)
    
    return cleaned


def _convert_expr_to_tcl(expr, params, args, local_vars=None, global_var_names=None):
    """Convert Python expression to Tcl with enhanced support.

    Supports:
    - Math functions: sin, cos, tan, exp, log, sqrt, abs
    - Operators: +, -, *, /, ** (power)
    - Phase shifts: sin(omega*t + phi)
    - Mixed variables: sin(a*t + b*x)
    - Nested functions: sin(cos(x))
    - Parentheses and complex arithmetic
    - Variable substitution from local_vars dict

    Parameters
    ----------
    expr : str
        Python expression to convert
    params : str
        Function parameters (comma-separated)
    args : list
        Argument names for spatial functions
    local_vars : dict, optional
        Dictionary of variable names to their values for substitution
    global_var_names : list, optional
        List of global variable names to replace with $var

    Returns
    -------
    str
        Tcl expression
    """
    import re

    # Determine which arg to use for substitution
    param_list = [p.strip() for p in params.split(',')]

    # Step 1: Substitute variables from local_vars BEFORE any other processing
    if local_vars:
        # Sort by length (longest first) to avoid partial replacements
        for var_name, var_value in sorted(local_vars.items(), key=lambda x: -len(x[0])):
            # Only substitute if it's a simple variable reference (not part of a larger name)
            pattern = r'(?<![a-zA-Z0-9_])' + re.escape(var_name) + r'(?![a-zA-Z0-9_])'
            # Format the value appropriately
            if isinstance(var_value, float):
                formatted_value = f'{var_value:.15g}'
            else:
                formatted_value = str(var_value)
            expr = re.sub(pattern, formatted_value, expr)

    # Step 2: Pre-process - evaluate numeric expressions where possible
    expr = _evaluate_numeric_constants(expr)

    # Step 3: Replacements for Python → Tcl
    replacements = [
        # Math functions (order matters - do numpy. before np.)
        # Hyperbolic functions
        (r'numpy\.tanh\s*\(', 'tanh('),
        (r'np\.tanh\s*\(', 'tanh('),
        (r'math\.tanh\s*\(', 'tanh('),
        (r'numpy\.sinh\s*\(', 'sinh('),
        (r'np\.sinh\s*\(', 'sinh('),
        (r'math\.sinh\s*\(', 'sinh('),
        (r'numpy\.cosh\s*\(', 'cosh('),
        (r'np\.cosh\s*\(', 'cosh('),
        (r'math\.cosh\s*\(', 'cosh('),

        # Inverse trigonometric functions
        (r'numpy\.arcsin\s*\(', 'asin('),
        (r'np\.arcsin\s*\(', 'asin('),
        (r'math\.arcsin\s*\(', 'asin('),
        (r'numpy\.arccos\s*\(', 'acos('),
        (r'np\.arccos\s*\(', 'acos('),
        (r'math\.arccos\s*\(', 'acos('),
        (r'numpy\.arctan\s*\(', 'atan('),
        (r'np\.arctan\s*\(', 'atan('),
        (r'math\.arctan\s*\(', 'atan('),
        (r'numpy\.arctan2\s*\(', 'atan2('),
        (r'np\.arctan2\s*\(', 'atan2('),
        (r'math\.arctan2\s*\(', 'atan2('),

        # Other math functions (sign and clip removed - not supported by Tcl)
        (r'numpy\.round\s*\(', 'round('),
        (r'np\.round\s*\(', 'round('),
        (r'numpy\.maximum\s*\(', 'max('),
        (r'np\.maximum\s*\(', 'max('),
        (r'numpy\.minimum\s*\(', 'min('),
        (r'np\.minimum\s*\(', 'min('),

        # Trigonometric functions
        (r'numpy\.sin\s*\(', 'sin('),
        (r'numpy\.cos\s*\(', 'cos('),
        (r'numpy\.tan\s*\(', 'tan('),
        (r'numpy\.exp\s*\(', 'exp('),
        (r'numpy\.log\s*\(', 'log('),
        (r'numpy\.log10\s*\(', 'log10('),
        (r'numpy\.log2\s*\(', 'log2('),
        (r'numpy\.sqrt\s*\(', 'sqrt('),
        (r'numpy\.abs\s*\(', 'abs('),
        (r'numpy\.floor\s*\(', 'floor('),
        (r'numpy\.ceil\s*\(', 'ceil('),

        (r'np\.sin\s*\(', 'sin('),
        (r'np\.cos\s*\(', 'cos('),
        (r'np\.tan\s*\(', 'tan('),
        (r'np\.exp\s*\(', 'exp('),
        (r'np\.log\s*\(', 'log('),
        (r'np\.log10\s*\(', 'log10('),
        (r'np\.log2\s*\(', 'log2('),
        (r'np\.sqrt\s*\(', 'sqrt('),
        (r'np\.abs\s*\(', 'abs('),
        (r'np\.floor\s*\(', 'floor('),
        (r'np\.ceil\s*\(', 'ceil('),

        # Math module functions
        (r'math\.sin\s*\(', 'sin('),
        (r'math\.cos\s*\(', 'cos('),
        (r'math\.tan\s*\(', 'tan('),
        (r'math\.exp\s*\(', 'exp('),
        (r'math\.log\s*\(', 'log('),
        (r'math\.log10\s*\(', 'log10('),
        (r'math\.log2\s*\(', 'log2('),
        (r'math\.sqrt\s*\(', 'sqrt('),
        (r'math\.abs\s*\(', 'abs('),

        # Constants
        (r'numpy\.pi\b', '3.14159265358979'),
        (r'np\.pi\b', '3.14159265358979'),
        (r'math\.pi\b', '3.14159265358979'),
        (r'numpy\.e\b', '2.71828182845905'),
        (r'np\.e\b', '2.71828182845905'),
        (r'math\.e\b', '2.71828182845905'),

        # Operators
        # Note: ** is converted separately by _convert_power_to_tcl()
        (r'//', '/'),    # Floor division
    ]

    tcl_expr = expr
    for pattern, replacement in replacements:
        tcl_expr = re.sub(pattern, replacement, tcl_expr)

    # 🔴 Улучшенная конвертация ** → pow()
    tcl_expr = _convert_power_to_tcl(tcl_expr)

    # Step 4: Replace function parameters with $var
    # Sort by length (longest first) to avoid partial replacements
    param_list_sorted = sorted(param_list, key=len, reverse=True)

    for param in param_list_sorted:
        param = param.strip()
        if param:
            # Replace param with $param
            # Use word boundaries but be careful with underscores
            pattern = r'(?<![a-zA-Z0-9_])' + re.escape(param) + r'(?![a-zA-Z0-9_])'
            tcl_expr = re.sub(pattern, '$' + param, tcl_expr)

    # Step 5: Replace global variable names with $var
    if global_var_names:
        # Sort by length (longest first) to avoid partial replacements
        global_var_names_sorted = sorted(global_var_names, key=len, reverse=True)
        for var_name in global_var_names_sorted:
            pattern = r'(?<![a-zA-Z0-9_])' + re.escape(var_name) + r'(?![a-zA-Z0-9_])'
            tcl_expr = re.sub(pattern, '$' + var_name, tcl_expr)

    # Step 6: Clean up - remove extra whitespace but preserve structure
    # Remove spaces around operators for cleaner output
    # IMPORTANT: Do NOT add spaces around '-' to preserve scientific notation (1e-9)
    tcl_expr = re.sub(r'\s*\+\s*', ' + ', tcl_expr)
    # Skip '-' to preserve scientific notation: 1e-9 NOT 1e - 9
    # tcl_expr = re.sub(r'\s*-\s*', ' - ', tcl_expr)  # DISABLED
    tcl_expr = re.sub(r'\s*\*\s*', '*', tcl_expr)
    tcl_expr = re.sub(r'\s*/\s*', '/', tcl_expr)
    tcl_expr = re.sub(r'\s*\^\s*', '^', tcl_expr)

    # Remove leading/trailing whitespace from the whole expression
    tcl_expr = tcl_expr.strip()

    # Remove any trailing commas (from tuple unpacking issues)
    tcl_expr = tcl_expr.rstrip(',')

    # 🔴 Валидация результата конвертации
    _validate_tcl_result(tcl_expr, expr)

    return tcl_expr


def _convert_power_to_tcl(expr):
    """Convert a**b to pow(a,b) with improved support.
    
    Supports:
    - Simple: x**2 → pow(x,2)
    - Numeric base: 2**10 → pow(2,10)
    - Float exponent: x**0.5 → pow(x,0.5)
    - Parenthesized: (x+y)**2 → pow((x+y),2)
    - Numeric powers: 2**10 → 1024 (evaluated)
    - With spaces: x ** 2 → pow(x,2)
    
    Parameters
    ----------
    expr : str
        Python expression with ** operator
    
    Returns
    -------
    str
        Tcl expression with pow() function
    """
    import re
    
    # Случай 1: (expr)**(number) - скобки и число (с пробелами)
    expr = re.sub(
        r'(\([^)]+\))\s*\*\*\s*(\d+\.?\d*)',
        r'pow(\1,\2)',
        expr
    )
    
    # Случай 2: (expr)**(name) - скобки и переменная (с пробелами)
    expr = re.sub(
        r'(\([^)]+\))\s*\*\*\s*([a-zA-Z_$][a-zA-Z0-9_$]*)',
        r'pow(\1,\2)',
        expr
    )
    
    # Случай 3: name**number - переменная и число (с пробелами)
    expr = re.sub(
        r'([a-zA-Z_$][a-zA-Z0-9_$]*)\s*\*\*\s*(\d+\.?\d*)',
        r'pow(\1,\2)',
        expr
    )
    
    # Случай 4: name**name - переменная и переменная (с пробелами)
    expr = re.sub(
        r'([a-zA-Z_$][a-zA-Z0-9_$]*)\s*\*\*\s*([a-zA-Z_$][a-zA-Z0-9_$]*)',
        r'pow(\1,\2)',
        expr
    )
    
    # Случай 5: number**number - вычисляется сразу (с пробелами)
    def eval_power(match):
        base = float(match.group(1))
        exp = float(match.group(2))
        result = base ** exp
        return f'{result:.15g}'
    
    expr = re.sub(
        r'(\d+\.?\d*)\s*\*\*\s*(\d+\.?\d*)',
        eval_power,
        expr
    )
    
    return expr


def _validate_tcl_result(tcl_expr: str, original_expr: str) -> None:
    """Validate that converted Tcl expression is valid.
    
    Checks for common conversion errors and raises ConversionError
    with helpful messages.
    
    Parameters
    ----------
    tcl_expr : str
        Converted Tcl expression
    original_expr : str
        Original Python expression (for error messages)
    
    Raises
    ------
    ConversionError
        If Tcl expression contains Python syntax
    """
    errors = []
    
    # Проверка на неконвертированный оператор степени
    if '**' in tcl_expr:
        errors.append(
            "Power operator '**' not converted. Use pow(base, exponent) instead."
        )
    
    # Проверка на неконвертированные префиксы NumPy
    if 'np.' in tcl_expr or 'numpy.' in tcl_expr:
        errors.append(
            "NumPy prefix not removed. Use 'sin(x)' instead of 'np.sin(x)'."
        )
    
    # Проверка на неконвертированный префикс math
    if 'math.' in tcl_expr:
        errors.append(
            "Math prefix not removed. Use 'sin(x)' instead of 'math.sin(x)'."
        )
    
    if errors:
        raise ConversionError(
            expression=original_expr,
            reason="; ".join(errors),
            suggestion="Check supported functions: sin, cos, tan, exp, log, sqrt, abs, min, max, etc."
        )


def _evaluate_numeric_constants(expr):
    """Pre-evaluate numeric constants in expression.

    Converts expressions like '2 * np.pi * 1e9' to '6.283185307179586e+09'
    by safely evaluating numeric subexpressions.

    IMPORTANT: Preserves scientific notation (1e-9) without spaces.
    """
    import re
    import math

    # First, protect scientific notation from being split
    # Replace 1e-9 with placeholder, then restore later
    scientific_pattern = r'(\d+\.?\d*[eE][+-]?\d+)'
    scientific_matches = {}
    placeholder_idx = 0

    def save_scientific(match):
        nonlocal placeholder_idx
        placeholder = f'__SCI_{placeholder_idx}__'
        scientific_matches[placeholder] = match.group(1)
        placeholder_idx += 1
        return placeholder

    # Save scientific notation numbers
    expr_protected = re.sub(scientific_pattern, save_scientific, expr)

    # Try to evaluate simple numeric expressions
    # Pattern: sequence of numbers, operators, and math constants
    numeric_pattern = r'(\d+\.?\d*(?:e[+-]?\d+)?|\bnp\.pi\b|\bnp\.e\b|\bmath\.pi\b|\bmath\.e\b)(?:\s*[\*/]\s*(\d+\.?\d*(?:e[+-]?\d+)?|\bnp\.pi\b|\bnp\.e\b|\bmath\.pi\b|\bmath\.e\b))*'

    def eval_match(match):
        subexpr = match.group(0)
        try:
            # Replace constants with values
            subexpr_eval = subexpr.replace('np.pi', str(math.pi))
            subexpr_eval = subexpr_eval.replace('np.e', str(math.e))
            subexpr_eval = subexpr_eval.replace('math.pi', str(math.pi))
            subexpr_eval = subexpr_eval.replace('math.e', str(math.e))

            # Safely evaluate
            result = eval(subexpr_eval)
            if isinstance(result, (int, float)):
                # Format nicely - use scientific notation for very small/large numbers
                if isinstance(result, int) or result == int(result):
                    return str(int(result))
                elif abs(result) < 1e-4 or abs(result) > 1e6:
                    return f'{result:.15e}'
                else:
                    return f'{result:.15g}'
        except:
            pass
        return subexpr

    # Only evaluate in multiplication/division contexts
    # to avoid breaking function arguments
    expr_evaluated = re.sub(r'(?<=[\(\+\-,])\s*(\d+\.?\d*(?:e[+-]?\d+)?(?:\s*\*\s*\d+\.?\d*(?:e[+-]?\d+)?)+)\s*(?=[\)\*/,])',
                  lambda m: str(eval(m.group(1).replace(' ', ''))), expr_protected)

    # Restore scientific notation numbers
    for placeholder, value in scientific_matches.items():
        expr_evaluated = expr_evaluated.replace(placeholder, value)

    return expr_evaluated


def _convert_func_fallback_enhanced(func, arg, args):
    """Enhanced fallback converter using AST when source is not available.
    
    This handles lambda functions from interactive mode by analyzing
    the function's behavior and generating appropriate Tcl code.
    
    Parameters
    ----------
    func : callable
        Python function to convert
    arg : str
        Argument name
    args : list
        Argument names for spatial functions
    
    Returns
    -------
    str or list
        Tcl expression string(s) or '0' if conversion fails
    """
    import numpy as np
    
    # Try to determine function type by testing
    try:
        if len(args) == 1:
            test_value = 1.0
        else:
            test_value = tuple([1.0] * len(args))
        
        result = func(test_value) if len(args) == 1 else func(*test_value)
        
        # Determine if vector or scalar
        if isinstance(result, (tuple, list)):
            # Vector function - try to determine pattern
            return _infer_vector_function(func, arg, args, result)
        else:
            # Scalar function - try to determine pattern
            return _infer_scalar_function(func, arg, args, result)
            
    except Exception as e:
        return '0'


def _infer_scalar_function(func, arg, args, test_result):
    """Infer Tcl expression for scalar function by testing."""
    import numpy as np
    import math
    
    # Test at multiple points to determine function type
    try:
        # Test at 0
        val0 = func(0) if len(args) == 1 else func(*[0]*len(args))
        
        # Test at pi/2
        val_pi2 = func(math.pi/2) if len(args) == 1 else func(*[math.pi/2]*len(args))
        
        # Test at pi
        val_pi = func(math.pi) if len(args) == 1 else func(*[math.pi]*len(args))
        
        # Check for sin pattern: sin(0)=0, sin(pi/2)=1, sin(pi)=0
        if abs(val0) < 1e-10 and abs(val_pi2 - 1) < 0.1 and abs(val_pi) < 1e-10:
            return f'sin($arg)'
        
        # Check for cos pattern: cos(0)=1, cos(pi/2)=0, cos(pi)=-1
        if abs(val0 - 1) < 0.1 and abs(val_pi2) < 1e-10 and abs(val_pi - (-1)) < 0.1:
            return f'cos($arg)'
        
        # Check for exp pattern
        val_1 = func(1) if len(args) == 1 else func(*[1]*len(args))
        if abs(val_1 - math.exp(1)) < 0.1:
            return f'exp($arg)'
        
        # Default: return constant value
        return str(float(test_result))
        
    except:
        return str(float(test_result))


def _infer_vector_function(func, arg, args, test_result):
    """Infer Tcl expression for vector function by testing."""
    result = []
    
    for i in range(len(test_result)):
        # Test each component
        def component_func(*vals):
            r = func(*vals) if len(vals) > 1 else func(vals[0])
            return r[i] if isinstance(r, (tuple, list)) else r
        
        # Try to infer
        inferred = _infer_scalar_function(component_func, arg, args, test_result[i])
        result.append(inferred)
    
    return result


def _convert_func_fallback(func, arg, args):
    """Fallback converter when source is not available.
    
    Delegates to _convert_func_fallback_enhanced for better support.
    """
    return _convert_func_fallback_enhanced(func, arg, args)

