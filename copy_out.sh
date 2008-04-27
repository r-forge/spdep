REPOS=`pwd`
TARGET=`echo ${REPOS} | sed -e 's/spdep/svn_build/'`
echo ${REPOS} ${TARGET}
cd "${REPOS}/pkg"
rsync -auv --exclude=".svn" * "${TARGET}/spdep"

