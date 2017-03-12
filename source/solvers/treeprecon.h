#ifndef __TREEPRECON_H__
#define __TREEPRECON_H__

#include "../util/common.h"
#include "../util/matrix.h"
#include "../util/graph.h"
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
