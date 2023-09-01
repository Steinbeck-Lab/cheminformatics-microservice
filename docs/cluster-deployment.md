---
outline: deep
---

# Kubernetes deployment (K8S)

## Helm Charts - v0.1.6

Cheminformatics Python Microservices comes packaged with a [Helm](https://helm.sh/docs/) chart, that makes it easy to deploy and manage (scale) containers on [Kubernetes](https://kubernetes.io/) via a convenient package manager interface. 

By following the steps outlined in this documentation, you can easily deploy this microservice using [Helm](https://helm.sh/docs/).

For more information about Helm Charts based deployment please refer to official [Helm documentation](https://helm.sh/docs/).

**Prerequisites:**

Before proceeding with the deployment, ensure that you have the following:

Kubernetes cluster: Set up a functioning Kubernetes cluster with kubectl configured to interact with the cluster.

Helm: [Install Helm](https://helm.sh/docs/docs/intro/install/) on your local machine or the machine from which you'll be deploying the application. 

## Deployment using CPM Helm Chart
*  **Add repo:** Once Helm has been set up correctly, add the repo as follows: 
```bash
    helm repo add repo-helm-charts https://nfdi4chem.github.io/repo-helm-charts/
```
If you had already added this repo earlier, run `helm repo update` to retrieve the latest versions of the packages. You can then run `helm search repo repo-helm-charts` to see the charts.

* **Deploy the Chart:** To deploy the chart, use the helm install command followed by the chart package name and an optional release name:
```bash
helm install myrelease repo-helm-charts/cheminfo-microservice
```
The release name (myrelease in this example) is used to identify the deployment, and it must be unique within the Kubernetes cluster.
The above command with deploy the service with the default configuration provided in [values.yml](https://github.com/NFDI4Chem/repo-helm-charts/blob/main/charts/cheminfo-microservice/values.yaml) file. To overwrite the default configuration please follow this [link](https://helm.sh/docs/chart_template_guide/values_files/) to learn more.

* Helm will install the chart and deploy the application to your Kubernetes cluster. You can view the deployed resources using kubectl commands:
```bash
kubectl get pods    
kubectl get services
```
* **Upgrading and Managing Deployments:** To upgrade an existing deployment, Use the `helm upgrade` command to apply the changes to the existing release e.g.
```bash
helm upgrade myrelease repo-helm-charts/cheminfo-microservice-0.1.6
```

* **Uninstalling the Chart:** To remove a deployed chart and associated resources, use the helm uninstall command:
```bash
helm uninstall myrelease
````
This will delete all resources created by the chart, including pods, services, and any other Kubernetes objects.



## Google Kubernetes Engine(GKE)

This documentation will guide you through the process of deploying the CPM application in [Google Kubernetes Engine (GKE)](https://cloud.google.com/kubernetes-engine) using [Helm](https://helm.sh/docs/). Here we provide a step-by-step instructions to help you set up your environment, install Helm, and dpeloy the CPM application in the Google Kubernetes Cluster.

### Prerequisites
Before you begin, ensure you have the following prerequisites in place:
* Google Cloud Platform (GCP) Account: You need a GCP account to create a GKE cluster and use other GCP services.
* In the Google Cloud console, on the project selector page, select or [create a Google Cloud project](https://cloud.google.com/resource-manager/docs/creating-managing-projects).
* Enable the Compute Engine, Artifact Registry, and Google Kubernetes Engine APIs.
* Activate Cloud Shell
    * Go to [Google Cloud Console](https://console.cloud.google.com/)
    * Click on the 'Activate Cloud Shell' button located at the top of the Google Cloud console window.

####  Step 1: Create a GKE Cluster
A GKE cluster consists of a pool of Compute Engine VM instances running [Kubernetes](https://kubernetes.io/), the open source cluster orchestration system that powers GKE.
1. Navigate to [Google Kubernetes Engine](https://console.cloud.google.com/kubernetes) page in Google Cloud console.
2. Click on the Create icon.
3. You can choose between `Autopilot` or `Standard` and click on `Configure`. But as the maintainence cost for Autopilot is expensive so we go for Standard option in this documentation. Click [here](https://cloud.google.com/kubernetes-engine/docs/resources/autopilot-standard-feature-comparison) to know the difference between Autopilot and Standard.
4. Provide the suitable name for your cluster and select the region from the dropdown list.
<p align="center">
  <img align="center" src="/docs/gke-1.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="70%">
</p>
5. Choose the compute engine machine type from the Nodes section under default-pool. For CPM the minimum machine type requirement will be e2-standard-4(4vCPU,16 GB memory). Also make sure you enable the cluster autoscaling if you are planning to scale your pod vertically.
<p align="center">
  <img align="center" src="/docs/gke-2.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="70%">
</p>

<p align="center">
  <img align="center" src="/docs/gke-3.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="70%">
</p>

6. Click on `Create` button and wait for the Cluster to be ready.

#### Step2: Connect to Cluster
1. After you cluster is created you can see the green check next to it. Once you see the check, click on the three dots next to your cluster to click on `Connect`.
2. Click on `Run in Cloud Shell` option and press enter.

#### Step3: Deploy via Helm Chart
1. Add the helm repo by running below command.
```bash
    helm repo add repo-helm-charts https://nfdi4chem.github.io/repo-helm-charts/
```
2. Run the below command to deploy the CPM app.
```bash
helm install myrelease repo-helm-charts/cheminfo-microservice
```
The release name (myrelease in this example) is used to identify the deployment, and it must be unique within the Kubernetes cluster.
The above command with deploy the service with the default configuration provided in [values.yml](https://github.com/NFDI4Chem/repo-helm-charts/blob/main/charts/cheminfo-microservice/values.yaml) file. To overwrite the default configuration please follow this [link](https://helm.sh/docs/chart_template_guide/values_files/) to learn more.

3. Helm will install the chart and deploy the application to your Kubernetes Cluster. To see the progress click on `Workloads` & `Services` tab. or run below commands in your Cloud Shell.
```bash
kubectl get pods    
kubectl get services
```
4. Once all services are deployed you can see green check next to each services. If not then check the logs to learn more.
5. To access your service you may have to expose it either by setting the type to `Load Balancer` under service in values.yml file or by an Ingress depending upon your requirement.
To learn more about how to configure Nginx Ingress you can click on the link [here](https://github.com/GoogleCloudPlatform/community/blob/master/archived/nginx-ingress-gke/index.md#deploying-the-nginx-ingress-controller-with-helm).

#### Step4: Clean up
To avoid incurring charges to your Google Cloud account for the resources used in this tutorial, either delete the project that contains the resources, or keep the project and delete resources by running below command.
```bash
helm delete myrelease
```

### Scaling
In case of performance issue, the CPM application can be scaled automatically using different strategies. The Helm chart for CPM comes packaged with Horizontal Pod Autoscaler. Horizontal Pod Autoscaler (HPA) is a Kubernetes feature that automatically adjusts the number of replica pods in a deployment or replica set based on observed CPU utilization or other select metrics. This allows your Kubernetes cluster to dynamically respond to changes in application load, ensuring optimal resource utilization and application availability.
You can set the `targetCPUUtilizationPercentage` and `targetMemoryUtilizationPercentage` values in values.yml file  according to your need and demand, which is the deciding factor to trigger the scaling. To learn more about different types of scaling in Kubernetes and GKE follow the official [documentation](https://cloud.google.com/kubernetes-engine/docs/concepts/cluster-autoscaler) of Google Cloud.

::: info

For Docker Compose based deployments follow the documentation here for [scaling](/docker.html#scaling).

::: 

### Contribute or Report an issue with CPM Helm Chart
Thank you for your valuable assistance in enhancing our deployment process. If you would like to contribute, kindly create a pull request in our [GitHub](https://github.com/NFDI4Chem/repo-helm-charts) repository. For any issues or bugs you have encountered, please feel free to create an [issue](https://github.com/NFDI4Chem/repo-helm-charts/issues) in the same or write to us at caffeine-devs@uni-jena.de. 
Your feedback is greatly appreciated.


References - https://cloud.google.com/docs
