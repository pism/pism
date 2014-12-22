%{
#include "IceGrid.hh"
%}

%extend pism::IceGrid
{
    %pythoncode {
    def points(self):
        """Iterate over tuples ``(i,j)`` of nodes owned by the current processor."""
        for i in xrange(self.xs(),self.xs()+self.xm()):
            for j in xrange(self.ys(),self.ys()+self.ym()):
                yield (i,j)
    def points_with_ghosts(self,nGhosts=0):
        for i in xrange(self.xs()-nGhosts,self.xs()+self.xm()+nGhosts):
            for j in xrange(self.ys()-nGhosts,self.ys()+self.ym()+nGhosts):
                yield (i,j)
    def coords(self):
        for i in xrange(self.xs(),self.xs()+self.xm()):
            for j in xrange(self.ys(),self.ys()+self.ym()):
                yield (i,j,self.x(i),self.y(j))
    }
}

/* This is needed to wrap IceGrid::get_dm() */
%shared_ptr(pism::PISMDM)
%ignore pism::PISMDM::operator DM;

%include "IceGrid.hh"
