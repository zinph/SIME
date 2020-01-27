# PKS Enumerator V2

HTML GUI for PKS enumerator V2 and flask.

## Setup
For better results, create a conda environment and activate it like:

```sh
conda create -n v2pks python
```

then:
```sh
conda activate v2pks
```

Then install all needed dependencies from `requirements.txt` in the following way:
```sh
pip install -r requirements.txt
```

## How to Use
Start Flask server :
```sh
python server.py
```

Then type the following link on the browser to load GUI.

[localhost:5000](http://localhost:5000)

Select the parameters according to what you need for your chemical library. Select Submit to generate.

The resulting files will found in LIBRARIES folder with the time stamp in the file names.
