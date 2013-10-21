SConscript(['lib/SConstruct', 'test/SConstruct', 'samples/SConstruct'])

env = Environment(TARFLAGS = '-cz --exclude ".*" --exclude "*.o" --exclude "*.tar.gz"')
tar_list = ['AUTHORS', 'COPYING', 'ChangeLog', 'LICENSE.txt', 'NEWS',
            'NEWS.jp.utf8.txt', 'README', 'README.jp.utf8.txt',
            'SConstruct', 'localSconsLib.py']
tar_list += Glob('include/MTToolBox/*.hpp')
tar_list += Glob('test/*.cpp')
tar_list += Glob('test/*.hpp')
tar_list += Glob('test/*.c')
tar_list += Glob('test/*.h')
tar_list += ['test/SConstruct']
tar_list += Glob('lib/*.a')
tar_list += Glob('lib/*.cpp')
tar_list += ['lib/SConstruct']
tar_list += Glob('docs/*')
tar_list += Glob('samples/*')
tar_list += Glob('doxygen/*.conf')
tar_list += Glob('doxygen/*.txt')
tar_list += Glob('samples/*')
tar_list += Glob('slide/*.tex')
env.Tar('MTToolbox.0.2.tar.gz', tar_list);


