Import ('env')

name = 'backbones-generator'
inc = env.Dir('.')
deps = ['prot-filer', 'mili', 'getoptpp']

src = env.Glob('src/*.cpp')
src.remove(env.File('src/petu.cpp'))
env.CreateObject('backbones-generator-objects', inc, src, deps)

deps += ['backbones-generator-objects']
env.CreateProgram(name, inc, ['src/petu.cpp'], deps)

test_combinations_env = env.Clone()
test_combinations_env.Append(CPPDEFINES={'COMBINATIONS_DEBUG' : '1'})
test_combinations_env.CreateSharedLibrary('bbgen-combinations-test', inc, [], src, deps)

