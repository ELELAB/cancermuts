how to run the Cancermuts Streamlit application:
1. source cancermuts_venv/bin/activate
2. run create_datebase.py
    - requires input data directroy with file organisation 
    datadir/protein_class/protein_name/run_type/date/
    - requires entry details file of the entries to be displayed 
    (protein_class, protein_name, run_type)
    - requires output destination folder of the metatables
3. run create_index.py
4. streamlit run main_page.py