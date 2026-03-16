# Stage 1: Build the React application
FROM node:24-alpine AS build

# Set working directory
WORKDIR /app

# Copy package.json, package-lock.json, and scripts (needed by postinstall)
COPY package*.json ./
COPY scripts/ ./scripts/

# Install dependencies (HUSKY=0 skips git hooks setup in Docker)
ENV HUSKY=0
RUN npm ci --silent

# Copy source code
COPY . .

# Set build arguments with defaults
ARG VITE_API_URL=http://localhost:8000/latest
ENV VITE_API_URL=$VITE_API_URL

# Build the application
RUN npm run build

# Stage 2: Serve the static files with Nginx
FROM nginx:alpine

# Copy the build output from the previous stage
COPY --from=build /app/dist /usr/share/nginx/html

# Create a custom nginx config for React Router support
RUN echo 'server { \
    listen 80; \
    server_name localhost; \
    root /usr/share/nginx/html; \
    index index.html; \
    location / { \
        try_files $uri $uri/ /index.html; \
    } \
    location ~* \.(js|css|png|jpg|jpeg|gif|ico|svg)$ { \
        expires 1y; \
        add_header Cache-Control "public, max-age=31536000, immutable"; \
    } \
}' > /etc/nginx/conf.d/default.conf

# Expose port 80
EXPOSE 80

# Start Nginx server
CMD ["nginx", "-g", "daemon off;"]
