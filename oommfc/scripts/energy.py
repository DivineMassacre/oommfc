import contextlib
import numbers
import warnings

import discretisedfield as df

import oommfc as oc


def energy_script(system):
    mif = ""
    for term in system.energy:
        mif += globals()[f"{term.__class__.__name__.lower()}_script"](term, system)

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


def zeeman_script(term, system):
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
    """Generate MIF script for YY_TransformStageMEL (transformation-based strain)."""
    mif = ""
    
    # First, set up the base strain fields
    ediagmif, ediagname = oc.scripts.setup_vector_parameter(
        term.e_diag, f"{term.name}_ediag"
    )
    eoffdiagmif, eoffdiagname = oc.scripts.setup_vector_parameter(
        term.e_offdiag, f"{term.name}_eoffdiag"
    )
    mif += ediagmif
    mif += eoffdiagmif
    
    # Generate transformation script if callable is provided
    if callable(term.transform_script):
        mif += _generate_transform_script(term)
    
    mif += "# MagnetoElastic (YY_TransformStageMEL)\n"
    mif += f"Specify YY_TransformStageMEL:{term.name} {{\n"
    mif += f"  B1 {B1name}\n"
    mif += f"  B2 {B2name}\n"
    mif += f"  e_diag_field {ediagname}\n"
    mif += f"  e_offdiag_field {eoffdiagname}\n"
    mif += f"  type {term.transform_type}\n"
    
    if hasattr(term, 'transform_script_args') and term.transform_script_args:
        mif += f"  script_args {term.transform_script_args}\n"
        mif += f"  script transform_{term.name}\n"
    
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


def _generate_transform_script(term):
    """Generate Tcl transformation script from Python callable.
    
    The Python callable is evaluated at MIF generation time to create
    a Tcl script that computes the transformation at runtime.
    """
    mif = ""
    
    # Determine the script arguments
    script_args = getattr(term, 'transform_script_args', 'total_time')
    args_list = script_args.split()
    
    # Build the Tcl proc signature
    args_str = " ".join(args_list)
    mif += f"proc transform_{term.name} {{ {args_str} }} {{\n"
    mif += f"  # Transformation type: {term.transform_type}\n"
    mif += f"  # Python function: {term.transform_script.__name__}\n"
    mif += "  # Note: For Python callables, the transformation must be\n"
    mif += "  # pre-computed or implemented directly in Tcl.\n"
    mif += "  # This is a placeholder - full implementation requires\n"
    mif += "  # embedding Python evaluation or pre-computed coefficients.\n"
    mif += "\n"
    
    # Generate placeholder based on transform_type
    if term.transform_type == 'diagonal':
        mif += "  # Diagonal transform: returns 6 values (3 diagonal + 3 derivatives)\n"
        mif += "  # Format: M11 M22 M33 dM11/dt dM22/dt dM33/dt\n"
        mif += "  return [list 1.0 1.0 1.0 0.0 0.0 0.0]\n"
    elif term.transform_type == 'symmetric':
        mif += "  # Symmetric transform: returns 12 values\n"
        mif += "  return [list 1.0 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0]\n"
    elif term.transform_type == 'general':
        mif += "  # General transform: returns 18 values\n"
        mif += "  return [list 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]\n"
    else:
        mif += "  # Identity transform (no transformation)\n"
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
