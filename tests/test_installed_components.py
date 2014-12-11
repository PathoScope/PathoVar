import unittest
import os
import subprocess

class TestCallExternalApplications(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        print("!?")

    def test_call_samtools_from_path(self):
        call = subprocess.Popen("samtools", stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        msg = call.communicate()
        self.assertTrue(call.returncode == 1, "Failed to call `samtools` from shell\n\n %r" % (msg,))

    def test_call_snpeff_from_environment(self):
        self.assertTrue("SNPEFF_PATH" in os.environ, "Could not find environmental variable $SNPEFF_PATH to find snpEff's install directory")
        snpeff_path = os.environ['SNPEFF_PATH'] 
        call = subprocess.Popen("java -Xmx2g -jar " + snpeff_path + os.sep + "snpEff.jar -h", stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE, shell=True)
        msg = call.communicate()
        self.assertTrue(call.returncode == 255, "Failed to call `%s` from shell using environmental variable $SNPEFF_PATH.\n\n %r" % 
            ("java -Xmx2g -jar " + snpeff_path + os.sep + "snpEff.jar", msg))

    def test_call_blastp_from_path(self):
        call = subprocess.Popen("blastp -h", stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        msg = call.communicate()
        self.assertTrue(call.returncode == 0, "Failed to call `blastp` from shell.\n\n %r" % (msg,))


if __name__ == '__main__':
    unittest.main()