%{
#include "util/IceGrid.hh"
%}

%extend pism::IceGrid
{
    %pythoncode {
    def points(self):
        """Iterate over tuples ``(i,j)`` of nodes owned by the current processor."""
        for i in range(self.xs(),self.xs()+self.xm()):
            for j in range(self.ys(),self.ys()+self.ym()):
                yield (i,j)
    def points_with_ghosts(self,nGhosts=0):
        for i in range(self.xs()-nGhosts,self.xs()+self.xm()+nGhosts):
            for j in range(self.ys()-nGhosts,self.ys()+self.ym()+nGhosts):
                yield (i,j)
    def coords(self):
        for i in range(self.xs(),self.xs()+self.xm()):
            for j in range(self.ys(),self.ys()+self.ym()):
                yield (i,j,self.x(i),self.y(j))
    }
}

%shared_ptr(pism::IceGrid);
%include "util/IceGrid.hh"
