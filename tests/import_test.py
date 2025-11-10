def test_import() -> None:
    try:
        import my_project  # noqa: PLC0415
    except ImportError:
        my_project = None

    assert my_project is not None, "my_project module should not be None"
