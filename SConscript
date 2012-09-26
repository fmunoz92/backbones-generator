Import ('env')

name = 'backbones-generator'
inc = env.Dir('.')
deps = ['prot-filer', 'mili', 'getoptpp']

src = env.Glob('src/*.cpp')
src.remove(env.File('src/petu.cpp'))
env.CreateObject('backbones-generator-objects', inc, src, deps)

deps += ['backbones-generator-objects']
env.CreateProgram(name, inc, ['src/petu.cpp'], deps)
