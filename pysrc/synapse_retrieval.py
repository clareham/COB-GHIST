import synapseclient, synapseutils
import os
import gzip
import shutil
import tarfile
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()
syn = synapseclient.Synapse()
syn.login(authToken=os.getenv("SYNAPSE-PAT"))

# pull data and deposit in the local data directory
data_dir = Path(__file__).resolve().parent.parent / "data"
data_dir.mkdir(parents=True, exist_ok=True)
files = synapseutils.syncFromSynapse(syn, 'syn74933315', path=str(data_dir))

for gz_file in data_dir.rglob("*.gz"):
	if gz_file.suffixes[-2:] == [".tar", ".gz"]:
		with tarfile.open(gz_file, "r:gz") as tar:
			tar.extractall(path=gz_file.parent)
	else:
		output_file = gz_file.with_suffix("")
		with gzip.open(gz_file, "rb") as source, output_file.open("wb") as destination:
			shutil.copyfileobj(source, destination)
