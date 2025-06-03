---
outline: deep
---

# Local development setup

## ✅ Prerequisites

Before you begin, install the following:

1. **Git** - [https://git-scm.com/downloads](https://git-scm.com/downloads)
2. **Docker Desktop** - [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop)

---

## 🚀 Step 1: Clone the Repository

1. Open a terminal (Command Prompt on Windows, Terminal on Mac/Linux).
2. Run the following command to download the project:
   ```bash
   git clone https://github.com/Steinbeck-Lab/cheminformatics-microservice.git
   ```
3. Move into the project folder:
   ```bash
   cd cheminformatics-microservice
   ```

---

## 🐳 Step 2: Run the Application with Docker

1. Make sure Docker is running on your machine.
2. Use the below docker-compose file to start the project
```yaml
services:
  api:
    build:
      context: ./
      dockerfile: Dockerfile
    container_name: cheminformatics-microservice-api
    environment:
      - HOMEPAGE_URL=https://docs.api.naturalproducts.net
      - RELEASE_VERSION=v2.6.0
      - WORKERS=2
    ports:
      - "8000:80"
    volumes:
      - ./app:/code/app
      - ./requirements.txt:/code/requirements.txt
    healthcheck:
      test: curl -f http://localhost:8000/latest/chem/health || exit 1
      interval: 90s
      timeout: 10s
      retries: 20
      start_period: 60s
    restart: unless-stopped
    security_opt:
      - no-new-privileges:true
    networks:
      - cm_network

  web:
    build:
      context: ./frontend
      dockerfile: ./Dockerfile
    container_name: cheminformatics-microservice-frontend
    environment:
      - REACT_APP_API_URL=http://localhost:8000/latest
      - NODE_ENV=development
    ports:
      - "3000:80"
    volumes:
      - ./frontend/src:/app/src
      - ./frontend/public:/app/public
      - ./frontend/package.json:/app/package.json
      - ./frontend/package-lock.json:/app/package-lock.json
    depends_on:
      - api
    restart: unless-stopped
    networks:
      - cm_network

  prometheus:
    image: prom/prometheus:latest
    container_name: prometheus
    ports:
      - "9090:9090"
    volumes:
      - ./prometheus_data/prometheus.yml:/etc/prometheus/prometheus.yml:ro
      - prometheus_data:/prometheus
    command:
      - '--config.file=/etc/prometheus/prometheus.yml'
      - '--storage.tsdb.path=/prometheus'
      - '--web.console.libraries=/usr/share/prometheus/console_libraries'
      - '--web.console.templates=/usr/share/prometheus/consoles'
    restart: unless-stopped
    networks:
      - cm_network

  grafana:
    image: grafana/grafana:latest
    container_name: grafana
    user: "472"
    ports:
      - "3001:3000"
    volumes:
      - grafana_data:/var/lib/grafana
    environment:
      - GF_SECURITY_ADMIN_PASSWORD=${GRAFANA_ADMIN_PASSWORD:-admin}
      - GF_USERS_ALLOW_SIGN_UP=false
    restart: unless-stopped
    depends_on:
      - prometheus
    networks:
      - cm_network

volumes:
  prometheus_data:
    driver: local
  grafana_data:
    driver: local

networks:
  cm_network:
    name: cm_network
    driver: bridge

```
2. Run the following command to start the project:
   ```bash
   docker-compose up
   ```
3. Wait for Docker to finish setting up the containers.

---

## 🌐 Step 3: Access the Services

- **Backend API**: [http://localhost:8000](http://localhost:8000)
- **Frontend App**: [http://localhost:3000](http://localhost:3000)
- **Grafana Dashboard** (Metrics): [http://localhost:3001](http://localhost:3001)

---

## 📊 Step 4: Set Up Grafana Dashboard (First Time Only)

1. Visit [http://localhost:3001](http://localhost:3001)
2. Login with:
   - Username: `admin`
   - Password: `admin`
3. Change the password when prompted.
4. Click **“Add your first data source”** and choose **Prometheus**.
5. Set the URL as:
   ```
   http://prometheus:9090
   ```
6. Click **“Save & Test”**.
7. Download and import the prebuilt dashboard JSON from:
   [cm-dashboard.json](https://github.com/Steinbeck-Lab/cheminformatics-microservice/blob/main/cm-dashboard.json)

---
You may update the docker-compose file to disable or add additional services but by default, the docker-compose file shipped with the project has the backend and the front-end, Prometheus and Grafana (logging and visualisation of metrics) services and associated volumes shared via a network.


## Workers

Uvicorn also has the option to start and run several worker processes.

Nevertheless, as of now, Uvicorn's capabilities for handling worker processes are more limited than Gunicorn's. So, if you want to have a process manager at this level (at the Python level), then it might be better to try Gunicorn as the process manager.

In any case, you would run it like this:

<div class="termy">

```console
$ uvicorn main:app --host 0.0.0.0 --port 8080 --workers 4
```

Update the Dockerfile to watch for code changes

```
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--workers", "4", "--reload"]
```

</div>


## Logging (Prometheus and  Grafana)

::: info

The following instructions are based on the blog post - https://dev.to/ken_mwaura1/getting-started-monitoring-a-fastapi-app-with-Grafana-and-Prometheus-a-step-by-step-guide-3fbn

To learn more about using Grafana in general, see the official [Prometheus](https://prometheus.io/docs/introduction/overview/) and [Grafana](https://grafana.com/docs/) documentation, or check out our other monitoring tutorials.

:::


Prometheus and Grafana are useful tools to monitor and visualize metrics in FastAPI applications.

Prometheus is a powerful monitoring system that collects and stores time-series data. By instrumenting your FastAPI app with Prometheus, you can collect various metrics such as request count, response time, error rate, and resource utilization. Grafana is a popular data visualization tool that integrates seamlessly with Prometheus. It allows you to create custom dashboards and visualize the collected metrics in a meaningful and interactive way. With Grafana, you can build visual representations of your FastAPI app's performance, monitor trends, and gain insights into your application's behaviour.

CM docker-compose file comes prepackaged with Prometheus and Grafana services for you. When you run the docker-compose file these services also spin up automatically and will be available for you to monitor your application performance.

When you install CM for the first time you need to configure your Prometheus source and enable it as the Grafana data source. You can then use the data source to create dashboards.

### Grafana Dashboard
Now that we have Prometheus running we can create a Grafana dashboard to visualize the metrics from our FastAPI app. To create a Grafana dashboard we need to do the following:

1. Create a new Grafana dashboard.
2. Add a new Prometheus data source.
3. Add a new graph panel.
4. Add a new query to the graph panel.
5. Apply the changes to the graph panel.
6. Save the dashboard.
7. View the dashboard.
8. Repeat steps 3-7 for each metric you want to visualize.
9. Repeat steps 2-8 for each dashboard you want to create.
10. Repeat steps 1-9 for each app you want to monitor.

Once you have Grafana running go to: localhost:3000. You should see the following screen:

Grafana login

Enter the default username and password (admin/admin) and click "Log In". You should be prompted to change the password. Enter a new password and click "Save". You should see the following screen:

<p align="center">
  <img align="center" src="/docs/grafana_login.jpeg" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="auto">
</p>

Grafana home

Click on the "Create your first data source" button. You should see the following screen:

<p align="center">
  <img align="center" src="/docs/grafana.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="auto">
</p>

Grafana added the data source

<p align="center">
  <img align="center" src="/docs/grafana_ds.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="auto">
</p>


Click on the "Prometheus" button. You should see the following screen:

<p align="center">
  <img align="center" src="/docs/prometheus.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="auto">
</p>

Enter the following information:

Name: Prometheus <br/>
URL: http://Prometheus:9090 <br/>
Access: Server (Default) <br/>
Scrape interval: 15s <br/>
HTTP Method: GET <br/>
HTTP Auth: None <br/>
Basic Auth: None <br/>
With Credentials: No <br/>
TLS Client Auth: None <br/>
TLS CA Certificate: None <br/>

Click on the "Save & Test" button. You should see the following screen:

<p align="center">
  <img align="center" src="/docs/grafana_ds_saved.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="auto">
</p>

Click on the "Dashboards" button. You should see the following screen:

<p align="center">
  <img align="center" src="/docs/grafana_db.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="auto">
</p>

Click on the ""New Dashboard" button. You should see the following screen:

<p align="center">
  <img align="center" src="/docs/grafana_db_new.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="auto">
</p>

Download the Cheminformatics Microservice dashboard template (JSON) here - https://github.com/Steinbeck-Lab/cheminformatics-microservice/blob/main/cm-dashboard.json

## Benchmarking / Stress testing

[Vegeta](https://github.com/tsenart/vegeta) is an open-source command-line tool written in Go, primarily used for load testing and benchmarking HTTP services. It allows you to simulate a high volume of requests to a target URL and measure the performance characteristics of the service under various loads.


To perform stress testing using Vegeta, you can follow these steps:

1. Install Vegeta: Start by installing Vegeta on your machine. You can download the latest release binary from the official GitHub repository (https://github.com/tsenart/vegeta) or use a package manager like Homebrew (macOS/Linux) or Chocolatey (Windows) for installation.

2. Define a target endpoint: Identify the specific FastAPI endpoint you want to stress test. Make sure you have the URL and any necessary authentication or headers required to access the endpoint.

3. Prepare a Vegeta target file: Create a text file, e.g., `target.txt`, and define the target URL using the Vegeta target format. For example:

   ```plaintext
   GET http://localhost:8000/my-endpoint
   ```

   Replace `http://localhost:8000/my-endpoint` with the actual URL of your FastAPI endpoint.

4. Create a Vegeta attack plan: In another text file, e.g., `attack.txt`, specify the rate and duration for the stress test using the Vegeta attack format. For example:

   ```plaintext
   rate: 100
   duration: 10s
   ```

   This example sets the request rate to 100 requests per second for a duration of 10 seconds. Adjust the values according to your requirements.

5. Run the Vegeta attack: Open a terminal or command prompt, navigate to the directory where the target and attack files are located, and execute the following command:

   ```bash
   vegeta attack -targets=target.txt -rate=attack.txt | vegeta report
   ```

   This command runs the Vegeta attack using the target and attack files, sends requests to the specified FastAPI endpoint, and generates a report with statistics.

6. Analyze the stress test results: Vegeta will output detailed metrics and performance statistics for the stress test. It includes data on request rate, latency, success rate, and more. Analyze these results to evaluate the performance and stability of your FastAPI application under stress.

By following these steps, you can perform stress testing on your CM FASTAPI application using Vegeta, generating load and analyzing the performance characteristics of your endpoints. This process helps identify potential bottlenecks and validate the scalability of your application.

## Linting / Formatting

We recommend using Flake8 and Black to perform linting and formatting in Python

1. Install Flake8 and Black: Start by installing both Flake8 and Black. You can install them using pip by running the following command:

   ```bash
   pip install flake8 black
   ```

2. Linting with Flake8: Flake8 is a popular Python linter that checks your code for style and potential errors. Run flake8 by executing the following command in your project directory:

   ```bash
   flake8 --per-file-ignores="__init__.py:F401" --ignore E402,E501,W503 $(git ls-files '*.py') .
   ```

   Flake8 will analyze your code and provide feedback on any style violations or issues found.

3. Formatting with Black: Black enforces consistent code formatting and automatically formats your code. To format your code with Black, run the following command in your project directory:

   ```bash
   black .
   ```

   or

   ```bash
    black $(git ls-files '*.py') .
   ```

   The `.` specifies the current directory. Black will recursively format all Python files within the directory and apply the necessary formatting changes.

   Note: Black uses a strict formatting style, so it's a good practice to make sure you have committed your changes to a version control system before running Black.
