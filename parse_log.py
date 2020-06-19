import os

with open("log.dat", "r") as f:
  lines = f.readlines()

to_parse = ""
energy = []
proba = []
for tmp in lines:
  if( "Accepted" in tmp ):
    to_parse= tmp.replace("\n","").replace(" ","").split("#")[1]

    if( "P(n+1)" in to_parse ):
      proba.append(to_parse.split("=")[1])
      energy.append(to_parse.split("-")[0])
    else:
      energy.append(to_parse)

with open("energy.dat", "w") as f:
  for i in range(0, len(energy)):
    f.write(str(i)+ " " + energy[i] + "\n")

with open("proba.dat", "w") as f:
  for i in range(0, len(proba)):
    f.write(str(i)+ " " + proba[i] + "\n")


