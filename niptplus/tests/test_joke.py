from unittest import TestCase
import niptplus

class TestJoke(TestCase):
    def test_is_string(self):
        s = niptplus.joke()
        self.assertTrue(isinstance(s, basestring))
