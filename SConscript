Import ('env')

name = 'backbones-generator'
inc = env.Dir('.')
src = env.Glob('src/*.cpp')
deps = ['prot-filer', 'mili', 'getoptpp']

env.AppendUnique(CPPFLAGS = ['-DMILI_NAMESPACE'])
env.CreateProgram(name, inc, src, deps)

