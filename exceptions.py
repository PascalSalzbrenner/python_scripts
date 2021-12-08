# file to contain useful custom exception definitions

# written by Pascal Salzbrenner, pts28@cam.ac.uk

class InputError(Exception):
    """Exception raised when there is something generally wrong with the input"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
