# 1. Use Bullseye to completely bypass the host seccomp/GPG issue
FROM debian:bullseye-slim

# 2. Set your Miniconda paths
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# 3. Standard update and install will now work perfectly
RUN apt-get -o Acquire::Check-Valid-Until=false update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*

# ... [Your other installation steps go here] ...

# 4. Use WORKDIR instead of 'cd' to clear the original CMD warning
WORKDIR /PitViper/PitViper

# 5. The optimized exec form for stable container shutdowns
CMD ["bash", "run_docker.sh"]
