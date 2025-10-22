# CNT Research Tools Python Version
A resposiory for the python version of CNT-approved research tools.  

# Get Started 

Download or clone the toolbox into a local folder:
```
git clone git@github.com:penn-cnt/CNT_research_tools.git
```


## Python Set-Up

Dependencies: 
* [anaconda](https://www.anaconda.com)

Create a conda environment:
```
conda env create -n ieegpy -f python/ieegpy.yml
```
```
conda env create -n ieegpy pip
pip install -r requirements.txt
```

## Login Configuration
Generate bin password file *_ieeglogin.bin through
1. Matlab toolbox
2. Python provided create_pwd_file function

Put password file into the folder
Change config.json file to specify user name and password file name.

## Resources
* [Project Tracker Google doc](https://docs.google.com/spreadsheets/d/12f-cCzB2J7W96jZzbJH7HKbWrivrUa2PKcRpQZsHXpM/edit?usp=sharing)
* [Background Notes Google doc](https://docs.google.com/document/d/17qalWwt5yb7NOVwob53GO_U_6H2GWmNjsN_aptjJoSw/edit?usp=sharing)
