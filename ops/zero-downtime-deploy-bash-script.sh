#!/bin/bash

# Define variables
COMPOSE_FILE="/mnt/data/cheminformatics-python-microservice/ops/docker-compose-dev.yml"
DOCKER_REPO_NAME="nfdi4chem/cheminformatics-microservice:dev-latest"
IMAGE_NAME="nfdi4chem/cheminformatics-microservice:dev-latest"
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
                echo "Container not healthy.. Check again.."
         else
                IS_CONTAINER_HEALTHY=0
                break
         fi
   done

  # Remove old containers and images
  if [ $IS_CONTAINER_HEALTHY == 0 ] ; then
  
        # Set the desired container name prefix
        CONTAINER_NAME_PREFIX="ops_web"

        # Retrieve the container IDs that match the prefix
        container_ids=$(docker ps -a --filter "name=^/${CONTAINER_NAME_PREFIX}" --format "{{.ID}}")

        # Sort the container IDs by creation date in ascending order
        sorted_container_ids=$(echo "$container_ids" | xargs docker inspect --format='{{.Created}} {{.ID}}' | sort | awk '{print $2}')

        # Get the oldest container ID
        oldest_container_id=$(echo "$sorted_container_ids" | head -n 1)

        # Check if any container IDs were found
        if [[ -z "$oldest_container_id" ]]; then
                echo "No containers found with the name prefix '${CONTAINER_NAME_PREFIX}'."
                exit 1
        fi

        # Delete the old container and unused images
        docker stop $oldest_container_id
        docker rm $oldest_container_id
        docker image prune -af
        echo "Deleted the oldest container with ID: ${oldest_container_id}"

  else
        echo "Couldnot complete the deployment as the container is unhealthy.."
        docker stop $NEW_CONTAINER_ID
        docker rm $NEW_CONTAINER_ID
  fi

else
       echo "Skipping deployment as no new image available.."
fi
