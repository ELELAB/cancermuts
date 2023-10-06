
How to run the Cancermuts Streamlit application:

You will need following packages to run the application 
	- os
	- pandas
	- streamlit
	- matplotlib.pyplot
	- upsetplot

To run with the example database:
	1. streamlit run main_page.py

To run with extended database:
	1. Run create_database.py to retrieve the cancermuts metatables of all proteins of interest
		- Requires an entry details file with the information on the proteins to include in the database
                  with the following columns (protein_class, protein_name, run_type, date, source)
		- Requires input data directory
		- The cancermuts metatable must found in the following subdirectory of each protein
                  /datadir/protein_class/protein_name/run_type/date/
		- Requires output destination folder for the metatables
	2. Change the variable database_dir to the name of your data directory in the main_page.py file, e.g. 
	   database_dir = get_database_dir('NAME_OF_DATADIR')
        3. streamlit run main_page.py

