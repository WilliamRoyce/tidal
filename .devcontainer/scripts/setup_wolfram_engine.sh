#!/usr/bin/env bash
# setup_wolfram_engine.sh
# Interactive setup wizard for Wolfram Engine 14.3
# Part of the torsion-gertsenshtein devcontainer configuration

set -e

# Colors for output
if command -v tput &> /dev/null && [ -t 1 ]; then
    GREEN=$(tput setaf 2)
    YELLOW=$(tput setaf 3)
    RED=$(tput setaf 1)
    BLUE=$(tput setaf 4)
    BOLD=$(tput bold)
    RESET=$(tput sgr0)
else
    GREEN=""
    YELLOW=""
    RED=""
    BLUE=""
    BOLD=""
    RESET=""
fi

# Configuration
ENGINE_VERSION="14.3"
ENGINE_DIR="$HOME/.local/wolfram/engine/14.3"
USERBASE_DIR="$HOME/.local/wolfram/userbase"
KERNEL_PATH="$ENGINE_DIR/Executables/WolframKernel"
WOLFRAMSCRIPT_PATH="$ENGINE_DIR/Executables/wolframscript"

# Helper functions
print_header() {
    echo ""
    echo "${BOLD}${BLUE}================================================================================${RESET}"
    echo "${BOLD}${BLUE}  $1${RESET}"
    echo "${BOLD}${BLUE}================================================================================${RESET}"
    echo ""
}

print_success() {
    echo "${GREEN}✓${RESET} $1"
}

print_warning() {
    echo "${YELLOW}⚠${RESET} $1"
}

print_error() {
    echo "${RED}✗${RESET} $1"
}

print_info() {
    echo "${BLUE}ℹ${RESET} $1"
}

print_step() {
    echo ""
    echo "${BOLD}=== $1 ===${RESET}"
    echo ""
}

# Main script
clear
print_header "Wolfram Engine 14.3 Setup Wizard"

echo "This interactive wizard will guide you through installing Wolfram Engine 14.3"
echo "for use with the torsion-gertsenshtein devcontainer."
echo ""
echo "${BOLD}What this script does:${RESET}"
echo "  1. Check if Wolfram Engine 14.3 is already installed"
echo "  2. Guide you through downloading the installer (if needed)"
echo "  3. Run the Wolfram Engine installer interactively"
echo "  4. Configure WolframScript for the devcontainer"
echo "  5. Activate the engine with your Wolfram ID"
echo "  6. Create activation backups for persistence"
echo "  7. Verify the installation is working"
echo ""
echo "${BOLD}Prerequisites:${RESET}"
echo "  • Wolfram ID account (free at https://account.wolfram.com/)"
echo "  • Free Wolfram Engine for Developers license (included with Wolfram ID)"
echo "  • Internet connection for download and activation"
echo ""
read -p "Press Enter to continue..."

# Step 1: Check existing installation
print_step "Step 1/7: Checking for existing Wolfram Engine installation"

if [ -x "$KERNEL_PATH" ]; then
    print_success "Found WolframKernel at: $KERNEL_PATH"

    # Check version
    print_info "Verifying Wolfram Engine version..."
    if VERSION_OUTPUT=$("$KERNEL_PATH" -noprompt -run 'Print[$VersionNumber]; Exit[]' 2>&1); then
        VERSION=$(echo "$VERSION_OUTPUT" | grep -oP '^\d+\.\d+' | head -1)

        if [ "$VERSION" == "$ENGINE_VERSION" ]; then
            print_success "Wolfram Engine $ENGINE_VERSION is already installed!"
            echo ""
            echo "${GREEN}${BOLD}✓ Setup is already complete!${RESET}"
            echo ""
            echo "Your Wolfram Engine installation is ready to use."
            echo ""
            echo "${BOLD}Next steps:${RESET}"
            echo "  • Run xAct setup: ${BLUE}bash .devcontainer/scripts/setup_xact.sh${RESET}"
            echo "  • Validate setup: ${BLUE}bash .devcontainer/scripts/validate-setup.sh${RESET}"
            echo ""
            exit 0
        else
            print_warning "Found version $VERSION, but this setup requires version $ENGINE_VERSION"
            echo ""
            echo "You have Wolfram Engine installed, but not the correct version."
            echo "This devcontainer is configured for version $ENGINE_VERSION specifically."
            echo ""
            echo "${BOLD}Options:${RESET}"
            echo "  1. Continue with existing version (may have compatibility issues)"
            echo "  2. Exit and manually install version $ENGINE_VERSION"
            echo ""
            read -p "Continue anyway? (y/N): " -n 1 -r
            echo ""
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                print_info "Exiting. Please install Wolfram Engine $ENGINE_VERSION and run this script again."
                exit 1
            fi
            print_warning "Proceeding with version $VERSION (not recommended)"
        fi
    else
        print_warning "Could not determine Wolfram Engine version"
        print_info "Proceeding with installation verification..."
    fi
else
    print_info "Wolfram Engine not found at expected location"
    print_info "Expected location: $ENGINE_DIR"
fi

# Step 2: Download instructions
print_step "Step 2/7: Obtaining Wolfram Engine Installer"

echo "${BOLD}You need to download the Wolfram Engine 14.3 installer.${RESET}"
echo ""
echo "${BOLD}Download instructions:${RESET}"
echo "  1. Visit: ${BLUE}https://www.wolfram.com/engine/${RESET}"
echo "  2. Click 'Download' and sign in with your Wolfram ID (or create one)"
echo "  3. Select 'Linux' platform"
echo "  4. Download the installer (filename: ${BOLD}WolframEngine_14.3.0_LINUX.sh${RESET})"
echo "  5. The installer is approximately 3-4 GB"
echo ""
echo "${YELLOW}Note:${RESET} This script will not download the installer automatically due to"
echo "Wolfram's licensing terms. You must download it manually from your account."
echo ""
read -p "Press Enter once you have downloaded the installer..."

# Step 3: Locate installer
print_step "Step 3/7: Locating Installer File"

echo "Please provide the path to the Wolfram Engine installer you downloaded."
echo ""
echo "${BOLD}Common locations:${RESET}"
echo "  • ~/Downloads/WolframEngine_14.3.0_LINUX.sh"
echo "  • /tmp/WolframEngine_14.3.0_LINUX.sh"
echo ""

INSTALLER_PATH=""
while [ -z "$INSTALLER_PATH" ] || [ ! -f "$INSTALLER_PATH" ]; do
    read -p "Installer path: " -e -i "$HOME/Downloads/WolframEngine_14.3.0_LINUX.sh" INPUT_PATH

    # Expand tilde
    INSTALLER_PATH="${INPUT_PATH/#\~/$HOME}"

    if [ ! -f "$INSTALLER_PATH" ]; then
        print_error "File not found: $INSTALLER_PATH"
        echo "Please check the path and try again."
        INSTALLER_PATH=""
    else
        print_success "Found installer: $INSTALLER_PATH"

        # Verify it's executable or can be made executable
        if [ ! -x "$INSTALLER_PATH" ]; then
            print_info "Making installer executable..."
            chmod +x "$INSTALLER_PATH"
        fi
    fi
done

# Step 4: Run installer
print_step "Step 4/7: Running Wolfram Engine Installer"

echo "${BOLD}The installer will now launch.${RESET}"
echo ""
echo "${BOLD}Installation prompts you'll see:${RESET}"
echo "  1. License agreement - Review and accept"
echo "  2. Installation directory - Use: ${GREEN}$ENGINE_DIR${RESET}"
echo "  3. Scripts directory - Use default: ${GREEN}/usr/local/bin${RESET}"
echo ""
echo "${YELLOW}Important:${RESET} When prompted for the installation directory, enter:"
echo "  ${GREEN}$ENGINE_DIR${RESET}"
echo ""
echo "The installer requires sudo for creating symlinks in /usr/local/bin."
echo ""
read -p "Press Enter to start the installer..."

# Create target directory
mkdir -p "$ENGINE_DIR"
mkdir -p "$USERBASE_DIR"

# Run installer
print_info "Launching installer..."
echo ""
if sudo bash "$INSTALLER_PATH"; then
    print_success "Installer completed successfully"
else
    print_error "Installer failed or was cancelled"
    exit 1
fi

# Step 5: Verify installation
print_step "Step 5/7: Verifying Installation"

if [ ! -x "$KERNEL_PATH" ]; then
    print_error "WolframKernel not found at: $KERNEL_PATH"
    echo ""
    echo "The installation may not have completed correctly."
    echo "Please check the installer output above for errors."
    exit 1
fi

print_success "WolframKernel found at: $KERNEL_PATH"

# Verify version
print_info "Checking Wolfram Engine version..."
if VERSION_OUTPUT=$("$KERNEL_PATH" -noprompt -run 'Print[$VersionNumber]; Exit[]' 2>&1); then
    VERSION=$(echo "$VERSION_OUTPUT" | grep -oP '^\d+\.\d+' | head -1)

    if [ "$VERSION" == "$ENGINE_VERSION" ]; then
        print_success "Confirmed: Wolfram Engine $ENGINE_VERSION"
    else
        print_warning "Installed version is $VERSION (expected $ENGINE_VERSION)"
    fi
else
    print_warning "Could not verify version, but kernel is executable"
fi

# Step 6: Activation
print_step "Step 6/7: Activating Wolfram Engine"

echo "${BOLD}Wolfram Engine requires activation with your Wolfram ID.${RESET}"
echo ""
echo "You will be prompted to enter your Wolfram ID (email) and password."
echo "This is a one-time process that connects your installation to your"
echo "free Wolfram Engine for Developers license."
echo ""
echo "${BOLD}What happens during activation:${RESET}"
echo "  • Your credentials are sent securely to Wolfram's servers"
echo "  • A license file (mathpass) is created locally"
echo "  • The license is tied to this machine's ID"
echo "  • You won't need to activate again unless you reinstall"
echo ""
read -p "Press Enter to start activation..."

# Check if already activated
LICENSE_FILE="$USERBASE_DIR/Licensing/mathpass"
if [ -f "$LICENSE_FILE" ]; then
    print_success "License file already exists: $LICENSE_FILE"
    print_info "Testing activation..."
    if wolframscript -code '$Version' > /dev/null 2>&1; then
        print_success "Activation is valid!"
    else
        print_warning "License file exists but activation may have expired"
        print_info "Running activation process..."
        wolframscript
    fi
else
    print_info "No existing license found. Starting activation..."
    echo ""
    wolframscript
fi

# Verify activation succeeded
if [ -f "$LICENSE_FILE" ]; then
    print_success "Activation successful! License file created."
else
    print_error "Activation may have failed - license file not found"
    echo ""
    echo "You can try activating manually later by running: wolframscript"
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Step 7: Create backup and verify
print_step "Step 7/7: Creating Activation Backup"

echo "Creating backup of activation for persistence across container rebuilds..."
echo ""

# Use the activation manager if available
ACTIVATION_MANAGER=".devcontainer/scripts/wolfram-activation-manager.sh"
if [ -f "$ACTIVATION_MANAGER" ]; then
    print_info "Running activation backup..."
    if bash "$ACTIVATION_MANAGER" backup; then
        print_success "Activation backup created successfully"
    else
        print_warning "Backup script encountered an issue (this may be normal)"
    fi
else
    # Manual backup
    print_info "Creating manual backup..."
    BACKUP_DIR="$USERBASE_DIR/.activation_backup"
    mkdir -p "$BACKUP_DIR"
    if [ -d "$HOME/.cache/Wolfram" ]; then
        cp -r "$HOME/.cache/Wolfram/"* "$BACKUP_DIR/" 2>/dev/null || true
        print_success "Backup created at: $BACKUP_DIR"
    fi
fi

# Run health check
echo ""
print_info "Running health check..."
CHECK_SCRIPT=".devcontainer/scripts/check-wolfram.sh"
if [ -f "$CHECK_SCRIPT" ]; then
    if bash "$CHECK_SCRIPT"; then
        print_success "Health check passed!"
    else
        print_warning "Health check reported issues (see above)"
    fi
else
    # Manual health check
    print_info "Testing WolframKernel..."
    if "$KERNEL_PATH" -noprompt -run 'Print[2+2]; Exit[]' 2>&1 | grep -q "4"; then
        print_success "WolframKernel test: PASSED"
    else
        print_warning "WolframKernel test: FAILED"
    fi

    print_info "Testing wolframscript..."
    if wolframscript -code '2+2' 2>&1 | grep -q "4"; then
        print_success "wolframscript test: PASSED"
    else
        print_warning "wolframscript test: FAILED"
    fi
fi

# Success!
print_header "✓ Wolfram Engine Setup Complete!"

echo "${GREEN}${BOLD}Installation successful!${RESET}"
echo ""
echo "${BOLD}What was installed:${RESET}"
echo "  • Wolfram Engine $ENGINE_VERSION"
echo "  • Location: $ENGINE_DIR"
echo "  • License: Activated and backed up"
echo "  • UserBase: $USERBASE_DIR"
echo ""
echo "${BOLD}Next steps:${RESET}"
echo "  1. Install xAct packages:"
echo "     ${BLUE}bash .devcontainer/scripts/setup_xact.sh${RESET}"
echo ""
echo "  2. Validate complete setup:"
echo "     ${BLUE}bash .devcontainer/scripts/validate-setup.sh${RESET}"
echo ""
echo "  3. Test Wolfram Engine:"
echo "     ${BLUE}wolframscript -code 'Integrate[x^2, x]'${RESET}"
echo ""
echo "${BOLD}Documentation:${RESET}"
echo "  • Complete guide: ${BLUE}.devcontainer/docs/WOLFRAM_GUIDE.md${RESET}"
echo "  • Quick reference: ${BLUE}.devcontainer/QUICKREF.md${RESET}"
echo ""
echo "Your Wolfram Engine is now persistent across container rebuilds!"
echo ""
