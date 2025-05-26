#!/bin/bash

set -e

COMPOSE_FILE="/mnt/data/cheminformatics-microservice/ops/docker-compose-dev.yml"
API_IMAGE="nfdi4chem/cheminformatics-microservice:api-dev"
WEB_IMAGE="nfdi4chem/cheminformatics-microservice:app-dev"
NEW_CONTAINER_ID=""

# === Functions ===

check_health() {
    HEALTH=$(docker inspect --format='{{json .State.Health.Status}}' "$NEW_CONTAINER_ID")
    if [[ "$HEALTH" == *"healthy"* ]]; then
        echo "✅ Container is healthy."
        return 0
    else
        echo "⏳ Container is unhealthy or still starting."
        return 1
    fi
}

deploy_service(){
    local service_name=$1   # e.g., "api"

    IS_CONTAINER_HEALTHY=1

    echo "🚀 Scaling up $service_name container..."
    docker compose -f "$COMPOSE_FILE" up -d --scale "$service_name=2" --no-recreate

    NEW_CONTAINER_ID=$(docker ps -q -l)
    echo "🔍 New container ID: $NEW_CONTAINER_ID"

    echo "⏳ Waiting for new container to pass health check (up to 10 retries)..."
    for i in {1..10}; do
        if check_health; then
            IS_CONTAINER_HEALTHY=0
            break
        else
            echo "Retry $i/10: Waiting 60s..."
            sleep 60
        fi
    done

    if [ "$IS_CONTAINER_HEALTHY" == 0 ]; then
        echo "🧼 Replacing old $service_name container(s)..."

        # Retrieve and sort containers with matching name prefix
        container_ids=$(docker ps -a --filter "name=$service_name" --format "{{.ID}}")
        sorted_container_ids=$(echo "$container_ids" | xargs docker inspect --format='{{.Created}} {{.ID}}' | sort | awk '{print $2}')
        oldest_container_id=$(echo "$sorted_container_ids" | head -n 1)

        if [[ -z "$oldest_container_id" ]]; then
            echo "❌ No containers found with name prefix: ${service_name}"
            exit 1
        fi

        docker stop "$oldest_container_id"
        docker rm "$oldest_container_id"
        docker image prune -af

        echo "✅ Deleted old container ID: $oldest_container_id"
    else
        echo "❌ Deployment aborted: new $service_name container is unhealthy."
        docker stop "$NEW_CONTAINER_ID"
        docker rm "$NEW_CONTAINER_ID"
    fi
}

# === Image Update Check ===

echo "🔍 Checking for updated images..."

deploy_api=false
deploy_web=false

if [ "$(docker pull "$API_IMAGE" | grep -c "Status: Image is up to date")" -eq 0 ]; then
    echo "📦 New API image available."
    deploy_api=true
fi

if [ "$(docker pull "$WEB_IMAGE" | grep -c "Status: Image is up to date")" -eq 0 ]; then
    echo "📦 New WEB image available."
    deploy_web=true
fi

if [[ "$deploy_api" = false && "$deploy_web" = false ]]; then
    echo "✅ No updates. Deployment not needed."
    exit 0
fi

# === Deployment ===

if [[ "$deploy_api" = true ]]; then
    deploy_service api
fi

if [[ "$deploy_web" = true ]]; then
    deploy_service web
fi

echo "🎉 Deployment completed successfully!"