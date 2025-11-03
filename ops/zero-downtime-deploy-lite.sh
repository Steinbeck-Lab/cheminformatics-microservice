#!/bin/bash

set -e

COMPOSE_FILE="$(dirname "$0")/docker-compose.lite.yml"
API_IMAGE="nfdi4chem/cheminformatics-microservice:latest-lite"
NEW_CONTAINER_ID=""

# === Error Handling ===
trap 'echo "❌ Script failed at line $LINENO"; exit 1' ERR

# === Basic Logging ===
LOG_FILE="/var/log/cheminformatics-deploy.log"

# Ensure log directory exists
mkdir -p "$(dirname "$LOG_FILE")"

# Validate compose file exists and is readable
if [[ ! -f "$COMPOSE_FILE" ]]; then
    echo "❌ Compose file not found: $COMPOSE_FILE"
    exit 1
fi

log() {
    local message="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "$message"
    echo "$message" >> "$LOG_FILE"
}

# === Functions ===

check_health() {
    HEALTH=$(docker inspect --format='{{json .State.Health.Status}}' "$NEW_CONTAINER_ID")
    if [[ "$HEALTH" == *"healthy"* ]]; then
        log "✅ Container is healthy."
        return 0
    else
        log "⏳ Container is unhealthy or still starting."
        return 1
    fi
}

deploy_service(){
    local service_name=$1   # e.g., "api"

    IS_CONTAINER_HEALTHY=1

    log "🚀 Scaling up $service_name container..."
    docker compose -f "$COMPOSE_FILE" up -d --scale "$service_name=2" --no-recreate

    NEW_CONTAINER_ID=$(docker ps -q -l)
    log "🔍 New container ID: $NEW_CONTAINER_ID"

    log "⏳ Waiting for new container to pass health check (up to 10 retries)..."
    for i in {1..10}; do
        if check_health; then
            IS_CONTAINER_HEALTHY=0
            break
        else
            log "Retry $i/10: Waiting 60s..."
            sleep 60
        fi
    done

    if [ "$IS_CONTAINER_HEALTHY" == 0 ]; then
        log "🧼 Replacing old $service_name container(s)..."

        # Retrieve and sort containers with matching name prefix
        container_ids=$(docker ps -a --filter "name=$service_name" --format "{{.ID}}")
        sorted_container_ids=$(echo "$container_ids" | xargs docker inspect --format='{{.Created}} {{.ID}}' | sort | awk '{print $2}')
        oldest_container_id=$(echo "$sorted_container_ids" | head -n 1)

        if [[ -z "$oldest_container_id" ]]; then
            log "❌ No containers found with name prefix: ${service_name}"
            exit 1
        fi

        docker stop "$oldest_container_id"
        docker rm "$oldest_container_id"
        docker image prune -af

        log "✅ Deleted old container ID: $oldest_container_id"
    else
        log "❌ Deployment aborted: new $service_name container is unhealthy."
        docker stop "$NEW_CONTAINER_ID"
        docker rm "$NEW_CONTAINER_ID"
    fi
}

# === Image Update Check ===

log "🔍 Checking for updated images..."

deploy_api=false

if [ "$(docker pull "$API_IMAGE" | grep -c "Status: Image is up to date")" -eq 0 ]; then
    log "📦 New API image available."
    deploy_api=true
else
    log "✅ API image is already up to date."
fi

if [[ "$deploy_api" = false ]]; then
    log "✅ No updates. Deployment not needed."
    exit 0
fi
# === Deployment ===

if [[ "$deploy_api" = true ]]; then
    deploy_service api
fi

log "🎉 Deployment completed successfully!"