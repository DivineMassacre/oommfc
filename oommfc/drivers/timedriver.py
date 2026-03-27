from .driver import Driver


class TimeDriver(Driver):
    """Time driver.

    Only attributes in ``_allowed_attributes`` can be defined. For details on
    possible values for individual attributes and their default values, please
    refer to ``Oxs_TimeDriver`` documentation
    (https://math.nist.gov/oommf/doc/userguide21a0/userguidexml/sec_oxsDrivers.html).

    Attributes
    ----------
    stage_iteration_limit : int, optional
        Maximum number of iterations allowed per stage. If this limit
        is reached, the stage is considered done, even if other stopping
        criteria are not met.

        .. note::

            For spatiotemporal Zeeman simulations, this parameter is
            automatically set to 1 to ensure correct synchronization
            between stages and iterations. Without this setting, the
            simulation may produce incorrect results due to time
            desynchronization in Oxs_StageZeeman.

    Examples
    --------
    1. Defining driver with a keyword argument.

    >>> import oommfc as oc
    ...
    >>> td = oc.TimeDriver(total_iteration_limit=5)

    2. Passing an argument which is not allowed.

    >>> import oommfc as oc
    ...
    >>> td = oc.TimeDriver(myarg=1)
    Traceback (most recent call last):
       ...
    AttributeError: ...

    3. Getting the list of allowed attributes.

    >>> import oommfc as oc
    ...
    >>> td = oc.TimeDriver()
    >>> td._allowed_attributes
    [...]

    4. Spatiotemporal Zeeman (automatically sets stage_iteration_limit=1):

    >>> import oommfc as oc
    ...
    >>> driver = oc.TimeDriver()
    >>> # driver.drive(system, t=1e-9, n=200)  # stage_iteration_limit=1 set automatically

    5. Explicit stage_iteration_limit:

    >>> import oommfc as oc
    ...
    >>> driver = oc.TimeDriver()
    >>> driver.stage_iteration_limit = 1  # 1 iteration per stage
    >>> # driver.drive(system, t=1e-9, n=200)

    """

    _allowed_attributes = [
        "evolver",
        "stopping_dm_dt",
        "stage_iteration_limit",
        "total_iteration_limit",
        "stage_count_check",
        "checkpoint_file",
        "checkpoint_interval",
        "checkpoint_disposal",
        "start_iteration",
        "start_stage",
        "start_stage_iteration",
        "start_stage_start_time",
        "start_stage_elapsed_time",
        "start_last_timestep",
        "normalize_aveM_output",
        "report_max_spin_angle",
        "report_wall_time",
    ]

    def _checkargs(self, kwargs):
        t, n = kwargs["t"], kwargs["n"]
        if t <= 0:
            msg = f"Cannot drive with {t=}."
            raise ValueError(msg)
        if not isinstance(n, int):
            msg = f"Cannot drive with {type(n)=}."
            raise ValueError(msg)
        if n <= 0:
            msg = f"Cannot drive with {n=}."
            raise ValueError(msg)

    def _check_system(self, system):
        """Checks the system has dynamics in it"""
        if len(system.dynamics) == 0:
            raise RuntimeError("System's dynamics is not defined")
        if len(system.energy) == 0:
            raise RuntimeError("System's energy is not defined")

    @property
    def _x(self):
        return "t"
