"""Tests for magneto-elastic energy script generation."""

import pytest

import micromagneticmodel as mm
import oommfc as oc


class TestMagnetoElasticScript:
    """Tests for magnetoelastic MIF script generation."""

    def test_fixedmel_script_generation(self):
        """Test YY_FixedMEL script generation for static strain."""
        system = mm.System(name="test")

        # Create a simple mesh and magnetization
        import discretisedfield as df
        region = df.Region(p1=(0, 0, 0), p2=(10e-9, 10e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(2, 2, 2))

        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = mm.MagnetoElastic(
            B1=1e7, B2=1e7, e_diag=(1e-3, 1e-3, 1e-3), e_offdiag=(0, 0, 0)
        )

        # Generate script
        from oommfc.scripts.energy import magnetoelastic_script
        mel_term = system.energy.magnetoelastic
        mif = magnetoelastic_script(mel_term, system)

        # Check script content
        assert "YY_FixedMEL" in mif
        assert "B1" in mif
        assert "B2" in mif
        assert "e_diag_field" in mif
        assert "e_offdiag_field" in mif

    def test_stagemel_script_generation(self):
        """Test YY_StageMEL script generation for stage-based strain."""
        system = mm.System(name="test")

        import discretisedfield as df
        region = df.Region(p1=(0, 0, 0), p2=(10e-9, 10e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(2, 2, 2))

        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = mm.MagnetoElastic.stage(
            B1=1e7,
            B2=1e7,
            e_diag_files=["strain_0.ovf", "strain_1.ovf"],
            e_offdiag_files=["off_0.ovf", "off_1.ovf"],
            stage_count=2,
        )

        from oommfc.scripts.energy import magnetoelastic_script
        mel_term = system.energy.magnetoelastic
        mif = magnetoelastic_script(mel_term, system)

        # Check script content
        assert "YY_StageMEL" in mif
        assert "e_diag_files" in mif
        assert "strain_0.ovf" in mif
        assert "e_offdiag_files" in mif
        assert "stage_count" in mif
        assert "stage_count 2" in mif

    def test_transformstagemel_script_generation(self):
        """Test YY_TransformStageMEL script generation."""
        system = mm.System(name="test")

        import discretisedfield as df
        region = df.Region(p1=(0, 0, 0), p2=(10e-9, 10e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(2, 2, 2))

        def dummy_transform(t):
            return [1, 1, 1, 0, 0, 0]

        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = mm.MagnetoElastic.transform(
            B1=1e7,
            B2=1e7,
            e_diag=(1, 0.3, 0.3),
            e_offdiag=(0, 0, 0),
            transform_script=dummy_transform,
            transform_type="diagonal",
        )

        from oommfc.scripts.energy import magnetoelastic_script
        mel_term = system.energy.magnetoelastic
        mif = magnetoelastic_script(mel_term, system)

        # Check script content - should NOT have explicit script_args
        # (uses default: stage stage_time total_time)
        assert "YY_TransformStageMEL" in mif
        assert "type diagonal" in mif
        assert "proc transform_" in mif
        assert "e_diag_script" in mif
        assert "e_offdiag_script" in mif
        assert "proc strain_diag_" in mif
        assert "proc strain_offdiag_" in mif
        # Verify 3-argument signature
        assert "stage stage_time total_time" in mif
        # script field MUST be present
        assert "script transform_" in mif
        # script_args should NOT be present (uses default)
        assert "script_args" not in mif

    def test_transform_all_types(self):
        """Test all transformation types generate valid scripts."""
        import discretisedfield as df
        region = df.Region(p1=(0, 0, 0), p2=(10e-9, 10e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(2, 2, 2))

        for transform_type in ["identity", "diagonal", "symmetric", "general"]:
            system = mm.System(name=f"test_{transform_type}")
            system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
            system.energy = mm.MagnetoElastic.transform(
                B1=1e7,
                B2=1e7,
                e_diag=(1, 0.3, 0.3),
                e_offdiag=(0, 0, 0),
                transform_script=lambda t: [1] * 18,
                transform_type=transform_type,
            )

            from oommfc.scripts.energy import magnetoelastic_script
            mel_term = system.energy.magnetoelastic
            mif = magnetoelastic_script(mel_term, system)

            assert f"type {transform_type}" in mif

    def test_mel_class_property(self):
        """Test _mel_class property returns correct OOMMF class."""
        # Static
        mel_static = mm.MagnetoElastic.static(
            B1=1e7, B2=1e7, e_diag=(1e-3, 1e-3, 1e-3), e_offdiag=(0, 0, 0)
        )
        assert mel_static._mel_class == "YY_FixedMEL"

        # Stage-based
        mel_stage = mm.MagnetoElastic.stage(
            B1=1e7, B2=1e7,
            e_diag_files=["a.ovf"],
            e_offdiag_files=["b.ovf"]
        )
        assert mel_stage._mel_class == "YY_StageMEL"

        # Transform-based
        mel_transform = mm.MagnetoElastic.transform(
            B1=1e7, B2=1e7,
            e_diag=(1, 0.3, 0.3),
            e_offdiag=(0, 0, 0),
            transform_script=lambda t: [1] * 6
        )
        assert mel_transform._mel_class == "YY_TransformStageMEL"

    def test_transform_direct_substitution_mode(self):
        """Test direct substitution mode: func without e_diag/e_offdiag."""
        system = mm.System(name="test_direct")

        import discretisedfield as df
        region = df.Region(p1=(0, 0, 0), p2=(10e-9, 10e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(2, 2, 2))

        def strain_func(t):
            # Returns full strain values
            return [1e-3, 1e-3, 1e-3, 0, 0, 0]

        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = mm.MagnetoElastic.transform(
            B1=1e7,
            B2=1e7,
            func=strain_func,  # No e_diag/e_offdiag
            dt=1e-13,
        )

        from oommfc.scripts.energy import magnetoelastic_script
        mel_term = system.energy.magnetoelastic
        mif = magnetoelastic_script(mel_term, system)

        # Check direct substitution mode markers
        assert "YY_TransformStageMEL" in mif
        assert "Direct substitution mode" in mif
        assert "proc transform_" in mif
        # Base strain should be zero
        assert "vector { 0 0 0 }" in mif

    def test_transform_matrix_mode(self):
        """Test matrix transformation mode: func + e_diag/e_offdiag."""
        system = mm.System(name="test_matrix")

        import discretisedfield as df
        region = df.Region(p1=(0, 0, 0), p2=(10e-9, 10e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(2, 2, 2))

        def transform_matrix(t):
            # Returns transformation matrix M(t)
            return [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = mm.MagnetoElastic.transform(
            B1=1e7,
            B2=1e7,
            e_diag=(1e-3, 1e-3, 1e-3),  # Base strain
            e_offdiag=(0, 0, 0),
            func=transform_matrix,  # Returns M(t)
            dt=1e-13,
            transform_type='diagonal',
        )

        from oommfc.scripts.energy import magnetoelastic_script
        mel_term = system.energy.magnetoelastic
        mif = magnetoelastic_script(mel_term, system)

        # Check matrix transformation mode markers
        assert "YY_TransformStageMEL" in mif
        assert "Matrix transformation mode" in mif
        assert "e_final = M(t)" in mif
        # Base strain should be user-provided values
        assert "1e-03" in mif or "0.001" in mif

    def test_transform_tcl_strings_mode(self):
        """Test tcl_strings mode for advanced control."""
        system = mm.System(name="test_tcl")

        import discretisedfield as df
        region = df.Region(p1=(0, 0, 0), p2=(10e-9, 10e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(2, 2, 2))

        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = mm.MagnetoElastic(
            B1=1e7,
            B2=1e7,
            e_diag=(0, 0, 0),
            e_offdiag=(0, 0, 0),
            tcl_strings={'script': 'return [list 1e-3 1e-3 1e-3 0 0 0]'},
        )

        from oommfc.scripts.energy import magnetoelastic_script
        mel_term = system.energy.magnetoelastic
        mif = magnetoelastic_script(mel_term, system)

        # Check tcl_strings mode
        assert "YY_TransformStageMEL" in mif or "YY_FixedMEL" in mif
        assert "Custom tcl script" in mif or "tcl_strings" in mif
