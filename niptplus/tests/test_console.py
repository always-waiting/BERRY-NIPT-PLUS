# -*- coding: UTF-8 -*-
import sys
from nose import with_setup
from unittest import TestCase
from cStringIO import StringIO
from contextlib import contextmanager
from niptplus.command_line import niptplus_commands

@contextmanager
def stdout_redirect(where):
    sys.stdout = where
    try:
        yield where
    finally:
        sys.stdout = sys.__stdout__


def setup():
    print "命令行检测开始"

def teardown(self):
    print "命令行检测完毕"


class TestConsole(TestCase):

    @with_setup(setup, teardown)
    def test_basic(self):
        print "testing niptplus_commands() print output    --->    ",
        #with stdout_redirect(StringIO()) as new_stdout:
        #    niptplus_commands()
        #new_stdout.seek(0)
        #content = new_stdout.read();
        #self.assertEqual(content,"这是NIPTplus命令组\n")
        tmp = sys.argv[:]
        print tmp
        sys.argv = sys.argv[0:1]
        niptplus_commands()
        sys.argv = tmp
        print "done."
