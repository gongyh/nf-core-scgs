class panToolsException(Exception):
    """Base exception for all panTools exceptions thrown in this project."""

    def __init__(self, message=""):
        Exception.__init__(self, message)


class panToolsExit(Exception):
    """Raised when panTools is to quietly exit."""

    def __init__(self, message=""):
        Exception.__init__(self, message)


class incorrectDBPath(panToolsException):
    def __init__(self, message=""):
        panToolsException.__init__(self, message)


class rankNotFound(panToolsException):
    def __init__(self, message=""):
        panToolsException.__init__(self, message)


class pangenomeNotFound(panToolsException):
    def __init__(self, message=""):
        panToolsException.__init__(self, message)


class pathNotFound(panToolsException):
    def __init__(self, message=""):
        panToolsException.__init__(self, message)


class gffNotFound(panToolsException):
    def __init__(self, message=""):
        panToolsException.__init__(self, message)
