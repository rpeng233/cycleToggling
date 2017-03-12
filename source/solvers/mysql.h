#ifndef __MYSQL_H__
#define __MYSQL_H__

#include <mysql++.h>
#include "common.h"
#include "mysql.config.h"

using namespace mysqlpp;

char mysql_buffer[4000000];
char *mysql_table = NULL;

namespace Mysql
{
	Connection conn(false);
	Query *query;
	int CONNECT_SUCCESS = 0;
	int mysql_matid;
	vector<string> queries;
	
	void connect_database()
	{
		try {
			conn.connect(mysql_dbname, mysql_dbhost, mysql_username, mysql_password);
			query = new Query(conn.query());
		} catch (BadQuery er) { // handle any connection or                      
			// query errors that may come up 
			cerr << "Error: " << er.what() << endl; 
			exit(-1);
		} catch (const BadConversion& er) {
			// Handle bad conversions
			cerr << "Conversion error: " << er.what() << endl <<
				"\tretrieved data size: " << er.retrieved <<
				", actual size: " << er.actual_size << endl;
			exit(-1);
		} catch (const Exception& er) {
			// Catch-all for any other MySQL++ exceptions
			cerr << "Error: " << er.what() << endl;
			exit(-1);
		}
		CONNECT_SUCCESS = 1;
	}

	void updaterow()
	{
		if (!CONNECT_SUCCESS)
		{
			printf("%s\n",mysql_buffer); return;
		}
		stringstream ss;
		ss<<"update "<<mysql_table<<" set "<<mysql_buffer<<" where matid = '"<<mysql_matid<<"';";
		queries.push_back(ss.str());
	}
	
	void insertrow()
	{
		if (!CONNECT_SUCCESS)
		{
			printf("%s\n",mysql_buffer); return;
		}
		try {
			(*query)<<"insert into "<<mysql_table<<" "<<mysql_buffer<<";";
			SimpleResult res=query->execute();
			mysql_matid = res.insert_id();
		} catch (BadQuery er) { // handle any connection or                      
			// query errors that may come up 
			cerr << "Error: " << er.what() << endl; 
			exit(-1);
		} catch (const BadConversion& er) {
			// Handle bad conversions
			cerr << "Conversion error: " << er.what() << endl <<
				"\tretrieved data size: " << er.retrieved <<
				", actual size: " << er.actual_size << endl;
			exit(-1);
		} catch (const Exception& er) {
			// Catch-all for any other MySQL++ exceptions
			cerr << "Error: " << er.what() << endl;
			exit(-1);
		}
	}
	
	void doquery()
	{
		try {
			rept(it,queries)
			{
				(*query)<<(*it);
				query->execute();
			}
		} catch (BadQuery er) { // handle any connection or                      
			// query errors that may come up 
			cerr << "Error: " << er.what() << endl; 
			exit(-1);
		} catch (const BadConversion& er) {
			// Handle bad conversions
			cerr << "Conversion error: " << er.what() << endl <<
				"\tretrieved data size: " << er.retrieved <<
				", actual size: " << er.actual_size << endl;
			exit(-1);
		} catch (const Exception& er) {
			// Catch-all for any other MySQL++ exceptions
			cerr << "Error: " << er.what() << endl;
			exit(-1);
		}
	}
}

#endif
