docker build -t anammox-pipeline .
docker run --rm `
  -v ${PWD}:/work `
  anammox-pipeline
