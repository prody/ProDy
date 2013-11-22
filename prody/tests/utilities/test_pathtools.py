from prody.tests import TestCase

from glob import glob
from os import remove
from os.path import join

from prody.tests.datafiles import TEMPDIR
from prody.utilities import gunzip, openFile


class TestGunzip(TestCase):

    def setUp(self):

        self.pref = join(TEMPDIR, 'compressed.txt')
        self.gzfn = self.pref + '.gz'
        self.text = ''.join(['some random text '] * 100)
        try:
            self.bytes = bytes(self.text, encoding='utf-8')
        except TypeError:
            self.bytes = self.text
        out = openFile(self.gzfn, 'wt')
        out.write(self.text)
        out.close()

    def testFile(self):

        fn = gunzip(self.gzfn)
        text = open(fn).read()
        self.assertEqual(text, self.text)

    def testBuffer(self):

        buff = open(self.gzfn, 'rb').read()
        text = gunzip(buff)
        self.assertEqual(text, self.bytes)

    def tearDown(self):

        for fn in glob(self.pref + '*'):
            remove(fn)
