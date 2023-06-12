# reelGeneShiny
Shiny App from the hackathon of reelGene scores and other data

## Running on CBSU VM
### TLJH / Traefik Background
ReelGene is currently running on the "BucklerHub" JupyterHub server hosted on a VM in CBSU/BioHPC. To make use of the resources we already have, during hackathon this server was used to also host reelGene.

The BucklerHub runs [The Littlest JupyterHub](https://tljh.jupyter.org), which uses [Traefik](https://traefik.io/traefik/) as a reverse proxy to direct users to their own Jupyter servers that are spun up.

### Traefik Routing
TLJH allows you to add additional routes to Traefik by adding a `.toml` file in `/opt/tljh/state/rules` ([docs here](https://tljh.jupyter.org/en/latest/topic/escape-hatch.html#extending-rules-toml)). We utilize this feature to direct requests for https://reelgene.maizegenetics.net to the Shiny Docker containers. The specific config file we are using is [`reelgene.toml`](./reelgene.toml) in this repository (less the hashed password).

In order to keep reelGene private before publication, we are using [HTTP "Basic" Authentication](https://developer.mozilla.org/en-US/docs/Web/HTTP/Authentication) to provide a rudimentary login window before displaying the application. For example, if we wanted to login with username `spaghetti` and password `tasty`, we use `htpasswd` to generate a hash to enter into our config:

```
$ htpasswd -nbB spaghetti "tasty"
spaghetti:$2y$05$bajwHRzDUuGQWSkFayXZruhrGCYRvSdC7zfnl2vRlwohy5gRyhsCK
```

Copy the output and modify the file like so:

```toml
[frontends.reelgene.auth.basic]
users = [
	"spaghetti:$2y$05$bajwHRzDUuGQWSkFayXZruhrGCYRvSdC7zfnl2vRlwohy5gRyhsCK"
]
```

The configuration makes Traefik aware of four "servers" at ports 8001 - 8004 to which it distributes traffic in a round-robin fashion. `htpasswd` uses [`bcrypt`](https://en.wikipedia.org/wiki/Bcrypt) to hash passwords by default. It includes salt, such that generating a hash for the same combination results in unique hashes each time. In theory, it should be safe to commit this hashed password (but let's not do that).

### Running reelGene via Docker Compose
The four servers that load is distributed to correspond to four identical containers on the machine. These containers are managed through a simple Docker Compose file, [`docker-compose.yml`](./docker-compose.yml). To bring up the four containers just do:

```
docker compose up -d
```

The `-d` flag is to detach the containers from our current shell, so they run "in the background."

The containers themselves are built from this repo's [`Dockerfile`](./Dockerfile) which is based upon [`rocker/shiny-verse`](https://github.com/rocker-org/rocker-versioned2) with the additional required packages from CRAN and Bioconductor installed.