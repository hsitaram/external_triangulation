SOURCEDIR=../source
OBJFILES=${SOURCEDIR}/edge.o ${SOURCEDIR}/common.o ${SOURCEDIR}/polygon.o\
	 ${SOURCEDIR}/triangle.o ${SOURCEDIR}/triangulation.o ${SOURCEDIR}/trielement.o\
	 ${SOURCEDIR}/laplace_solver.o
OFILES=edge.o common.o polygon.o triangle.o triangulation.o trielement.o\
      laplace_solver.o
CXXFLAGS=-g -W -pedantic -Wall

dlny:	${OBJFILES}
	mv ${SOURCEDIR}/*.o .
	${CXX} ${OFILES} ${SOURCEDIR}/tri_driver.cpp -o dlny.exe
clean:
	rm -rf *.o *.exe *.dat
