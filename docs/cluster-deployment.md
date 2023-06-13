---
outline: deep
---

# Kubernetes deployment via Helm.
Cheminformatics Python Microservices comes packaged with [Helm](https://helm.sh/docs/) chart which simplifies the deployment and management of applications on [Kubernetes](https://kubernetes.io/) by providing a convenient package manager interface. By following the steps outlined in this documentation, you can easily deploy this microservice using [Helm](https://helm.sh/docs/), enabling efficient and reproducible deployments in your [Kubernetes](https://kubernetes.io/) cluster.

Please refer to Helm’s [documentation](https://helm.sh/docs/) to get started.

**Prerequisites:**

Before proceeding with the deployment, ensure that you have the following:

Kubernetes cluster: Set up a functioning Kubernetes cluster with kubectl configured to interact with the cluster.
Helm: [Install Helm](https://helm.sh/docs/docs/intro/install/) on your local machine or the machine from which you'll be deploying the application. 

### Deploy the chart:
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
helm upgrade myrelease repo-helm-charts/cheminfo-microservice-0.0.2
```

* **Uninstalling the Chart:** To remove a deployed chart and associated resources, use the helm uninstall command:
```bash
helm uninstall myrelease
````
This will delete all resources created by the chart, including pods, services, and any other Kubernetes objects.


### Contribute or Report an issue
Thank you for your valuable assistance in enhancing our deployment process. If you would like to contribute, kindly create a pull request in our [GitHub](https://github.com/NFDI4Chem/repo-helm-charts) repository. For any issues or bugs you have encountered, please feel free to create an [issue](https://github.com/NFDI4Chem/repo-helm-charts/issues) in the same or write to us at caffeine-devs@uni-jena.de. 
Your feedback is greatly appreciated.



