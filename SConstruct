# This is a comment


import os
import subprocess



AddOption('--cxx',
	dest='cxx',
	type='string',
	nargs=1,
	action='store',
	help='cxx compiler')


AddOption('--gprof',
	action="store_true",
	dest='gprof',
	default=False)


include = [Dir('/opt/local/include/'), '#/src/external/include/', '/usr/include/eigen3/', '/opt/local/include/eigen3/', '/home/jmacdona/dev/ceres-solver/include/']
lib_path = [Dir('/opt/local/lib/'), '/home/jmacdona/dev/ceres-solver/ceres-bin/lib/'];

libs = ['boost_system-mt', 'boost_filesystem-mt', 'boost_thread-mt', 'boost_timer-mt', 'libboost_chrono-mt', 'pthread', 'ceres', 'lapack', 'blas', 'gfortran']


cppdefines = ['USING_SCONS_BUILDER', ('HAVE_INLINE', 1), 'NDEBUG', ('BOOST_FILESYSTEM_VERSION', 3), 'BOOST_UBLAS_NDEBUG' ] 

gcc_cppflags = ['-msse3', '-funroll-loops', '-pipe','-g', '-Wall', '-fmessage-length=0', '-std=c++11' ]


# optimisation level normally -O3
gcc_cppflags = gcc_cppflags + ['-O3']


# FOR profiling
if GetOption('gprof') == True:
	gcc_cppflags = gcc_cppflags + ['-pg']
	print "INFO: compiling with gprof flags"


utils_src = Glob('src/utils/*.cpp')


common_sources = utils_src

gcc_link_flags = []
gcc_link_flags = gcc_link_flags + ['-static']
gcc_link_flags = gcc_link_flags + ['-g']
gcc_link_flags = gcc_link_flags + ['-rdynamic']


if GetOption('gprof') == True:
	gcc_link_flags = gcc_link_flags + ['-pg']
	print "INFO: linking with gprof flags"



gcc_env = Environment(ENV = os.environ, CC = 'gcc', CXX = 'g++', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags, LINKFLAGS=gcc_link_flags )   # Create an environmnet
gcc_env.Append(CPPPATH=["#/src"] + [include])

icc_env = Environment(ENV = os.environ, CC = 'icc', CXX = 'icpc', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags + ['-ip'], LINKFLAGS=gcc_link_flags + ['-ip'] )
icc_env.Append(CPPPATH=["#/src"] + [include])

clang_env = Environment(ENV = os.environ, CC = 'clang', CXX = 'clang++', LIBPATH = lib_path, LIBS = libs, CPPDEFINES = cppdefines, CPPFLAGS = gcc_cppflags, LINKFLAGS=gcc_link_flags )   # Create an environmnet
clang_env.Append(CPPPATH=["#/src"] + [include])


if GetOption('cxx') == None:
	print "--cxx option not set: ", GetOption('cxx')
	env = gcc_env
	print "setting env to gcc_env"
elif GetOption('cxx') == "g++" or GetOption('cxx') == "gcc":
        print "--cxx option set: ", GetOption('cxx')
        env = gcc_env
        print "setting env to gcc_env"
elif GetOption('cxx') == "clang++" or GetOption('cxx') == "clang":
	print "--cxx option set: ", GetOption('cxx')
	env = clang_env
	print "setting env to clang_env"
elif GetOption('cxx') == "icc" or GetOption('cxx') == "icpc":
        print "--cxx option set: ", GetOption('cxx')
        env = icc_env
        print "setting env to icc_env"
else:
        print "--cxx option set: ", GetOption('cxx')
	print "ERROR: C++ compiler option not recognised"



print "CC is:", env['CC']
print "CXX is:", env['CXX']
print "LINK is:", env['LINK']
print "TOOLS is:", env['TOOLS']

env.Program(target = "bin/ptmcmc_bmeg_xylose", source = ["src/apps/ptmcmc_bmeg_xylose.cpp"] + common_sources )
env.Program(target = "bin/ptmcmc_bmeg_competition", source = ["src/apps/ptmcmc_bmeg_competition.cpp"] + common_sources )
env.Program(target = "bin/mcmc_bmeg_xylose", source = ["src/apps/mcmc_bmeg_xylose.cpp"] + common_sources )
env.Program(target = "bin/mcmc_bmeg_competition", source = ["src/apps/mcmc_bmeg_competition.cpp"] + common_sources )



env.Program(target = "bin/output_trajectory_cond_def_competition", source = ["src/apps/output_trajectory_cond_def_competition.cpp"] + common_sources )
env.Program(target = "bin/output_trajectory_list_cond_def_competition", source = ["src/apps/output_trajectory_list_cond_def_competition.cpp"] + common_sources )
env.Program(target = "bin/output_trajectory_cond_def_xylose", source = ["src/apps/output_trajectory_cond_def_xylose.cpp"] + common_sources )
env.Program(target = "bin/output_trajectory_list_cond_def_xylose", source = ["src/apps/output_trajectory_list_cond_def_xylose.cpp"] + common_sources )
env.Program(target = "bin/output_trajectory_cond_def_competition_track", source = ["src/apps/output_trajectory_cond_def_competition_track.cpp"] + common_sources )
env.Program(target = "bin/output_trajectory_list_cond_def_competition_track", source = ["src/apps/output_trajectory_list_cond_def_competition_track.cpp"] + common_sources )
env.Program(target = "bin/output_2d_endpoint_matrix_competition", source = ["src/apps/output_2d_endpoint_matrix_competition.cpp"] + common_sources )
env.Program(target = "bin/output_2d_endpoint_matrix_xylose", source = ["src/apps/output_2d_endpoint_matrix_xylose.cpp"] + common_sources )




