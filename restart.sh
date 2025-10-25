# 1. Clear caches and restart services
sudo systemctl restart shiny-server
sudo systemctl restart nginx

# 2. Clear R temporary files
sudo rm -rf /tmp/Rtmp*
sudo rm -rf ~/.Rcache

# 3. Touch the app.R file to force reload
sudo touch /home/smith6jt/panc_CODEX/Islet-Explorer/app/shiny_app/app.R

# 4. If using symlinks, refresh them
# (adjust paths as needed for your /srv/shiny-server setup)




Development Mode (Local Testing)
Starting the App for Development

# Navigate to project directory
cd /home/smith6jt/panc_CODEX/Islet-Explorer

# Method 1: Direct R command
R -q -e "shiny::runApp('app/shiny_app', launch.browser=FALSE)"

Method 2: Using VS Code task (if available)
Run the "Run Shiny app to reproduce issues" task

Method 3: Interactive R session
R
> library(shiny)
> runApp('app/shiny_app', port=3838, host='127.0.0.1')


Deployment Mode (Shiny-Server + Nginx)
Current Deployment Architecture

Internet → Nginx (Port 8080) → Shiny-Server → App
          Reverse Proxy        Process Manager

Starting Deployment Services
# Start shiny-server (manages R processes)
sudo systemctl start shiny-server
sudo systemctl enable shiny-server  # Auto-start on boot

# Start nginx (reverse proxy for port 8080)
sudo systemctl start nginx
sudo systemctl enable nginx

# Check status
sudo systemctl status shiny-server
sudo systemctl status nginx

Stopping Deployment Services

# Stop services
sudo systemctl stop shiny-server
sudo systemctl stop nginx

# Or restart to pick up changes
sudo systemctl restart shiny-server
sudo systemctl restart nginx


When Code Changes Are Made:

# 1. Update the timestamp to force reload
sudo touch /home/smith6jt/panc_CODEX/Islet-Explorer/app/shiny_app/app.R

# 2. Clear caches (if persistent issues)
sudo rm -rf /tmp/Rtmp*
sudo rm -rf ~/.Rcache
sudo rm -rf /var/lib/shiny-server/bookmarks/*

# 3. Restart shiny-server to clear R process cache
sudo systemctl restart shiny-server

# 4. Optional: Clear nginx cache if needed
sudo rm -rf /var/cache/nginx/* 2>/dev/null || true
sudo systemctl restart nginx

Using the Automated Script:

# Use the cache clearing script we created
./scripts/clear_shiny_cache.sh

# Or run it with sudo for full service restart
sudo ./scripts/Shell_scripts/clear_shiny_cache.sh


App Access URLs
Development Mode

# Local development
http://127.0.0.1:3838/islet-explorer/
http://localhost:3838/islet-explorer/
# (port may vary - check terminal output)

Deployment Mode
# Remote access via nginx reverse proxy
http://10.15.152.7:8080/islet-explorer/
# The nginx configuration routes this to shiny-server internally

