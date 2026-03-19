"""Tests for spatiotemporal Zeeman MIF script generation."""

import numpy as np
import pytest

import micromagneticmodel as mm
import oommfc as oc
import discretisedfield as df


class TestZeemanSpatiotemporalScript:
    """Tests for spatiotemporal Zeeman MIF generation."""

    def test_uniform_time_term(self):
        """Test uniform time-dependent term."""
        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        zeeman = mm.Zeeman(H=(1e6, 0, 0))
        zeeman.add_time_term(lambda t: (1e3 * np.sin(2 * np.pi * 1e9 * t), 0, 0))

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = zeeman

        from oommfc.scripts.energy import zeeman_script
        mif = zeeman_script(zeeman, system)

        assert "FieldPerStage" in mif
        assert "SpatiotemporalField" in mif
        assert "Oxs_StageZeeman" in mif
        assert "current_time" in mif

    def test_gaussian_mask(self):
        """Test term with Gaussian spatial mask."""
        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        zeeman = mm.Zeeman()
        zeeman.add_time_term(
            func=lambda t: np.sin(2 * np.pi * 1e9 * t),
            mask=lambda x, y, z: (1e3 * np.exp(-x**2 / 1e-16), 0, 0)
        )

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = zeeman

        from oommfc.scripts.energy import zeeman_script
        mif = zeeman_script(zeeman, system)

        assert "exp" in mif
        assert "Oxs_StageZeeman" in mif

    def test_traveling_wave(self):
        """Test traveling wave (two terms)."""
        H0 = 1e6

        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        zeeman = mm.Zeeman()
        # Use simpler functions that the converter can handle
        zeeman.add_time_term(
            func=lambda t: H0 * np.sin(1e9 * t),
            mask=lambda x, y, z: np.cos(1e7 * x)
        )
        zeeman.add_time_term(
            func=lambda t: H0 * np.cos(1e9 * t),
            mask=lambda x, y, z: np.sin(1e7 * x)
        )

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = zeeman

        from oommfc.scripts.energy import zeeman_script
        mif = zeeman_script(zeeman, system)

        # Check that both terms are present
        assert mif.count("# Term") == 2

    def test_static_field_only(self):
        """Test that static field only uses FixedZeeman."""
        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        zeeman = mm.Zeeman(H=(1e6, 0, 0))

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = zeeman

        from oommfc.scripts.energy import zeeman_script
        mif = zeeman_script(zeeman, system)

        assert "Oxs_FixedZeeman" in mif
        assert "FieldPerStage" not in mif

    def test_combined_static_and_dynamic(self):
        """Test combining static field with dynamic term."""
        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        zeeman = mm.Zeeman(H=(1e6, 0, 0))
        zeeman.add_time_term(lambda t: (1e5 * np.sin(2 * np.pi * 1e9 * t), 0, 0))

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = zeeman

        from oommfc.scripts.energy import zeeman_script
        mif = zeeman_script(zeeman, system)

        # Check static field
        assert "H_static_x 1000000.0" in mif
        # Check dynamic term
        assert "FieldPerStage" in mif
        assert "Oxs_StageZeeman" in mif
