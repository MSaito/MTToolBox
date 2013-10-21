#
# Local Scons Library
#

#
# define check for libUnitTest++.a
#
def checkHeader(header, path):
    import os
    for a in path:
        if os.path.exists(os.path.join(a, header)):
            return True
    return False

def getSubPath(name, path, dirname):
    import os
    for f in path:
        if os.path.exists(os.path.join(f, name)):
            return f
        if os.path.exists(os.path.join(f, dirname, name)):
            return os.path.join(f, dirname)
    return False
#
# define unitTest
#
def runUnitTest(env, target, source):
    import subprocess
    app = str(source[0].abspath)
    if not subprocess.call(app):
        open(str(target[0]),'w').write("PASSED\n")

def getOSLibPath():
    import os
    if os.name == 'posix' or os.name == 'mac':
        os_library_path = ['/usr/local/lib', '/usr/local/lib64',
                           '/usr/lib', 'usr/lib64']
    else:
        os_library_path = []
    if os.environ.get('LIBRARY_PATH'):
        os_library_path += os.environ.get('LIBRARY_PATH').split(':')
    if os_library_path :
        while '' in os_library_path: os_library_path.remove('')
    path = []
    for f in os_library_path:
        if os.path.exists(f) and not (f in path) : path.append(f)
    return path

def getOSIncPath():
    import os
    if os.name == 'posix' or os.name == 'mac':
        os_include_path = ['/include', '/usr/include', '/usr/local/include']
    else:
        os_include_path = []
    if os.environ.get('CPLUS_INCLUDE_PATH'):
        os_include_path += os.environ.get('CPLUS_INCLUDE_PATH').split(':')
    path = []
    if os_include_path :
        while '' in os_include_path: os_include_path.remove('')
        for f in os_include_path:
            if os.path.exists(f) and not (f in path):
                path.append(f)
    return path

def convIncPath(in_path):
    import os
    path = []
    for f in in_path:
        element = '-I' + f
        if not element in path:
            path.append(element)
    return path
