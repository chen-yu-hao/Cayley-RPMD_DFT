import sys
if sys.version.startswith('3.1'):
    import unittest
else:
    import unittest2 as  unittest
suite = unittest.TestLoader().discover('tests', pattern='*test.py', top_level_dir='.')
unittest.TextTestRunner(verbosity=2).run(suite)
