#ifndef __MATLAB_H__
#define __MATLAB_H__

#include "common.h"

namespace Matlab
{
	string gettmpfile()
	{
		string file="/tmp/";
		rep(i,1,50) file+=('0'+rand()%10);
		return file;
	}
	
	void execute(string cmd)
	{
		system("rm /tmp/matlab.finish.lock > /dev/null 2&>1");
		cmd = "diary('/dev/null'); diary on; "+cmd+ " diary off; diary('/tmp/matlab.finish.lock'); diary on; disp (1); diary off;";
		ofstream fout("/tmp/matlab.pipe");
		fout<<cmd<<endl;
		fout.close();
		while (1)
		{
			FILE *f = fopen("/tmp/matlab.finish.lock","r");
			if (f)
			{
				fclose(f); break;
			}
		}
		system("rm /tmp/matlab.finish.lock");
	}
}

#endif 
