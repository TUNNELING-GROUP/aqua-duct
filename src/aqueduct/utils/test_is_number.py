from unittest import TestCase

def is_number(s):
    #http://pythoncentral.org/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
    if isinstance(s,bool):
        return False
    try:
        float(s)
        return True
    except:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except:
        pass
    return False

class TestIs_number(TestCase):

    def test_is_number(self):
        s= {1,"a","%"}
        Oczekiwany= {True, False, False}
        current =is_number(s)

        if ( current != Oczekiwany):
            Exception("test failed")
