#ifndef __TREEPRECON_H__
#define __TREEPRECON_H__

#include "common.h"
#include "matrix.h"
#include "graph.h"
#include "kosz.h"

function<Vec(const Vec&)> TreePreconditioner(const GraphSP &g)
{
	AbstractSolver S=KOSZ(g);
	return [S](const Vec &b) {
		Vec x; int flag; tie(x,flag)=S.solve(b,-1,0);
		return x;
	};
}

#endif
