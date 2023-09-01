---
outline: deep
---

# Sustainability and Maintainence

It's important to note that the cheminformatics toolkits distributed with CPM are all packaged under one container. We consciously chose to go against the usual notion/practice that containers are supposed to do one thing, so every cheminformatics toolkit needs to be packaged as a separate microservice. This is to avoid unnecessary complexity of container orchestration across multiple containers while the containers, as such, can scale indefinitely as they are stateless. 

This architectural preference might lead to compatiblity issues or pose a risk of conflicting package versions if the tools integrated are not actively maintained or upgraded (if the underlying toolkits have API changes).

While the microservice approach ensures that there will always be a working version of the package or tool accessible for anyone, ensuring reproducibility.

Incase of toolkit versioning compatibility issues we recommend switching those incompatible tools to dedicated docker image(s) and allow cross service interaction.

<hr/>

### Implementation

Install docker on the cheminformatics microservice (add the following line to Dockerfile)

```
RUN apt-get update && apt-get -y install docker.io
```

To let the docker client inside your CPM container interact with the tool image (docker service) on your host, you need to add /var/run/docker.sock as a volume to your CPM container. With Docker compose this is done by adding this to your docker service "volumes" section:

```
- /var/run/docker.sock:/var/run/docker.sock
```

Now when you build and run your docker images, you'll be able to execute "docker exec" from within the CPM docker, with a command like this, and you'll be talking to the tool image (docker service) on the host:

```
docker exec -u root <tool_service> /path/your_shell_script
```

From with in the CPM's FAST API Server you can use subprocess to execute a command on a different service and capture the output.

