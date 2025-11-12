TAG=ckhrulev/pism-intel

run: image
	docker run --rm -it ${TAG}

image: Dockerfile
	docker build -t ${TAG} --target pism-container -f $^ .

image-%: Dockerfile
	docker build -t ${TAG}-$* --target pism-$* -f $^ .

run-%: image-%
	docker run --rm -it ${TAG}-$*
