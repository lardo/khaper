# export
export PGDATA=/var/postgresql/data
export PERL=/lustre/project/og03/Public/Pipe/Software/PERL/perl_5.18.4/

# boost
export BOOST_LIB="/nfs/config/boost/lib"
export BOOST_INCLUDE="/nfs/config/boost/include"

# gcc
export GCC_LIB="/nfs2/config/gcc/gcc-4.9.2/lib/:/nfs2/config/gcc/gcc-4.9.2/lib64/"
export GXX="/nfs2/config/gcc/gcc-4.9.2/bin"
export CPLUS_INCLUDE_PATH="$BOOST_INCLUDE"
export MPC="/nfs2/lib/mpc-0.8.1/lib/"
export SEQ="/lustre/project/og03/Public/Git/TheThinker3/Source/DNA2"

# libarary
export LIB_BLASR="/lustre/project/og03/Public/Pipe/DNA_DENOVO/PacBio/pbdagcon/blasr_libcpp/alignment"
export LIBPBDATA_LIB="/lustre/project/og03/Public/Pipe/DNA_DENOVO/PacBio/pbdagcon/blasr_libcpp/pbdata"
export LD_LIBRARY_PATH="$LIB_BLASR:LIBPBDATA_LIB:$GCC_LIB:$BOOST_LIB:$MPC:$SEQ:$LD_LIBRARY_PATH"

# path
export PBDAGCON="/lustre/project/og03/Public/Pipe/DNA_DENOVO/PacBio/pbdagcon/src/cpp"
export PATH="$PBDAGCON:$GXX:$BOOST_INCLUDE:$PERL:$PATH"
