docker build -t anammox-pipeline .
docker run -it --rm `
  -v ${PWD}/reads:/work/reads `
  -v ${PWD}/hmm:/work/hmm `
  -v ${PWD}/results:/work/results `
  anammox-pipeline