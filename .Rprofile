# Ensure /usr/local/bin is in PATH for Pandoc and other tools
if (!grepl("/usr/local/bin", Sys.getenv("PATH"))) {
  Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep = ":"))
}

if (interactive()) {
  suppressMessages(require(devtools))
} 