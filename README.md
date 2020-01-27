# SIME

HTML GUI for SIME and flask.

## Setup
For better results, create a conda environment and activate it like:

```sh
conda create -n SIME python
```

then:
```sh
conda activate SIME
```

Then install all needed dependencies from `requirements.txt` in the following way:
```sh
pip install -r requirements.txt
```

## How to Use
Start Flask server :
```sh
python main.py
```

Then type the following line in your browser to load SIME software GUI.

[localhost:5000](http://localhost:5000)

Use the parameters as needed, and click **Submit** to generate your *in-silico* library of fully assembled macrolides.

The resulting info file and smile file(s) will be found in LIBRARIES folder with the time stamp within the file names. 
