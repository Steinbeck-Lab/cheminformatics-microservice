#!/bin/bash

# Define variables
COMPOSE_FILE="/mnt/data/cheminformatics-python-microservice/docker-compose-dev.yml"
DOCKER_REPO_NAME="caffeinejena/cheminformatics-python-microservice:dev-latest"
IMAGE_NAME="caffeinejena/cheminformatics-python-microservice:dev-latest"
NEW_CONTAINER_ID=""
IS_CONTAINER_HEALTHY=1

# Function to check the health of the container
check_health() {

    HEALTH=$(docker inspect --format='{{json .State.Health.Status}}' $NEW_CONTAINER_ID)

    if [[ $HEALTH == *"healthy"* ]]; then
        echo "Container is healthy."
        return 0
    else
        echo "Container is unhealthy or still starting"
        return 1
    fi

}

# Check if there is a new image available in the Docker repository
if [ "$(docker pull $DOCKER_REPO_NAME | grep "Status: Image is up to date" | wc -l)" -eq 0 ]; then
  # Pull the latest image
   echo "Pull the latest image.."
   docker-compose -f $COMPOSE_FILE pull

  # Scale up a new container
   echo "Scale up new container.."
   docker-compose -f $COMPOSE_FILE up -d --scale web=2 --no-recreate

   NEW_CONTAINER_ID=$(docker ps -q -l)

  echo "New Container Id is.."
  echo "$NEW_CONTAINER_ID"

  # Wait for new containers to start and health checks to pass
   echo "Waiting for the new containers to start and health check to pass retry 5 times.."
   n=0;
   while [ $n -le 10 ]
   do
         if ! check_health; then
                n=$(( $n + 1 ))
                sleep 1m
                echo "Container not healthy.. Check once more.."
         else
                IS_CONTAINER_HEALTHY=0
                break
         fi
   done

  # Remove old containers and images
   if [ $IS_CONTAINER_HEALTHY == 0 ] ; then

        # Get the container IDs of running containers with the same image
        container_ids=$(docker ps --filter ancestor=$IMAGE_NAME --format "{{.ID}}")

        # Check if there are exactly two containers running with the same image
        if [[ $(echo "$container_ids" | wc -l) -eq 2 ]]; then

                # Get the creation timestamps of the two containers
                timestamp1=$(docker inspect --format="{{.Created}}" $container_ids | head -n 1)
                timestamp2=$(docker inspect --format="{{.Created}}" $container_ids | tail -n 1)

                # Compare the timestamps to determine the older container
                if [[ "$timestamp1" < "$timestamp2" ]]; then
                        container_to_delete=$container_ids
                else
                        container_to_delete=$(echo "$container_ids" | tail -n 1)
                fi

                # Delete the older container and images
                docker stop $container_to_delete
                docker rm $container_to_delete
                docker image prune -af

                echo "Deleted container $container_to_delete"

                else
                        echo "Two running containers with the same image not found"
                fi

   else
        echo "Couldnot complete the deployment as the container is unhealthy.."
        docker stop $NEW_CONTAINER_ID
        docker rm $NEW_CONTAINER_ID
   fi

else
        echo "Skipping deployment as no new image available.."
fi