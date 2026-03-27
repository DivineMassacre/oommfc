"""Tests for stage_iteration_limit with spatiotemporal Zeeman."""

import numpy as np
import pytest
import warnings

import micromagneticmodel as mm
import oommfc as oc
import discretisedfield as df


class TestStageIterationLimit:
    """Tests for stage_iteration_limit=1 with spatiotemporal Zeeman."""

    def test_stage_iteration_limit_auto_set_for_spatiotemporal(self):
        """Test that stage_iteration_limit=1 is auto-set for spatiotemporal Zeeman."""
        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        zeeman = mm.Zeeman(H=(1e6, 0, 0))
        zeeman.add_time_term(lambda t: (1e3 * np.sin(2 * np.pi * 1e9 * t), 0, 0))

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = zeeman

        # Generate driver script
        driver = oc.TimeDriver()
        
        # Emulate driver.drive() behavior
        has_spatiotemporal = any(
            getattr(term, 'has_time_terms', False)
            for term in system.energy
        )
        if has_spatiotemporal and not hasattr(driver, 'stage_iteration_limit'):
            driver.stage_iteration_limit = 1
        
        mif = oc.scripts.driver_script(driver, system, t=1e-9, n=100)

        # Check stage_iteration_limit is set
        assert "stage_iteration_limit 1" in mif

    def test_stage_iteration_limit_warning(self):
        """Test that warning is issued for spatiotemporal without explicit stage_iteration_limit."""
        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        zeeman = mm.Zeeman(H=(1e6, 0, 0))
        zeeman.add_time_term(lambda t: (1e3 * np.sin(2 * np.pi * 1e9 * t), 0, 0))

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = zeeman

        driver = oc.TimeDriver()
        
        # Check warning is issued
        with pytest.warns(UserWarning, match="stage_iteration_limit=1"):
            has_spatiotemporal = any(
                getattr(term, 'has_time_terms', False)
                for term in system.energy
            )
            if has_spatiotemporal and not hasattr(driver, 'stage_iteration_limit'):
                warnings.warn(
                    "Spatiotemporal Zeeman requires stage_iteration_limit=1 "
                    "for correct time synchronization between stages and iterations. "
                    "Setting automatically. To disable this warning, explicitly set "
                    "driver.stage_iteration_limit = 1.",
                    UserWarning,
                    stacklevel=2
                )
                driver.stage_iteration_limit = 1

    def test_stage_iteration_limit_explicit(self):
        """Test that explicit stage_iteration_limit is preserved."""
        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = mm.Zeeman(H=(1e6, 0, 0))

        # Generate driver script with explicit stage_iteration_limit
        driver = oc.TimeDriver()
        driver.stage_iteration_limit = 10
        mif = oc.scripts.driver_script(driver, system, t=1e-9, n=100)

        # Check explicit value is used
        assert "stage_iteration_limit 10" in mif

    def test_no_stage_iteration_limit_for_static_zeeman(self):
        """Test that stage_iteration_limit is NOT set for static Zeeman."""
        region = df.Region(p1=(0, 0, 0), p2=(100e-9, 100e-9, 10e-9))
        mesh = df.Mesh(region=region, n=(10, 10, 2))

        system = mm.System(name="test")
        system.m = df.Field(mesh, nvdim=3, value=(1, 0, 0), norm=1)
        system.energy = mm.Zeeman(H=(1e6, 0, 0))

        # Generate driver script
        driver = oc.TimeDriver()
        mif = oc.scripts.driver_script(driver, system, t=1e-9, n=100)

        # stage_iteration_limit should NOT be in MIF for static Zeeman
        assert "stage_iteration_limit" not in mif
