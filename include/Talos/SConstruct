import distutils.sysconfig
env = Environment(SWIGFLAGS = ['-c++', '-python'],
                  CPPPATH = [distutils.sysconfig.get_python_inc()],
                  SHLIBPREFIX = "")

if ARGUMENTS["intel"] == "yes":
    print("Intel compiler is used in Talos.")
# if icpc is not found by Talos, please add Path to icpc.
#    env.AppendENVPath('PATH', '/ccc/products/icc-20.0.0/system/default/20.0.0/bin/intel64')
    env['CXX'] = 'icpc'

env.SharedLibrary('_talos.so', ['Talos.cpp', 'talos.i'])
