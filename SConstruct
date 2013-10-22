SConscript(['lib/SConstruct'])
version = '0.2'
tarfile = 'MTToolBox.' + version + '.tar.gz'
tardir = 'MTToolBox.' + version
tmpdir = '/tmp/MTToolBox.' + version
tarcommand = 'tar -cvzf ' + tarfile + ' -C/tmp '
tarcommand += '--exclude ".*" --exclude "*.o" --exclude "*.tar.gz" '
tarcommand += tardir
Command(tarfile, "SConstruct",
        [
        Copy(tmpdir, './'),
        tarcommand,
        Delete(tmpdir)
        ])
Default('lib')

