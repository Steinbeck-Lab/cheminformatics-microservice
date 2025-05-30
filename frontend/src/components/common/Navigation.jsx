import React from 'react';
import { NavLink } from 'react-router-dom';
import {
  HiOutlineHome,
  HiOutlineBeaker,
  HiOutlineRefresh, 
  HiOutlinePresentationChartLine,
  HiOutlineAdjustments,
  HiOutlineCamera,
  HiUserGroup
} from 'react-icons/hi';

// Add CSS for the strobing animation
const StrobingIcon = ({ icon: Icon, className }) => {
  return (
    <div className="relative">
      <Icon className={`${className} animate-usergroup-strobe`} aria-hidden="true" />
    </div>
  );
};

// Navigation links configuration
const navLinks = [
  {
    to: '/home',
    label: 'Home',
    icon: HiOutlineHome,
    exact: true,
  },
  {
    to: '/chem',
    label: 'Chemical Analysis',
    icon: HiOutlineBeaker,
  },
  {
    to: '/convert',
    label: 'Format Conversion',
    icon: HiOutlineRefresh,
  },
  {
    to: '/depict',
    label: 'Depiction',
    icon: HiOutlinePresentationChartLine,
  },
  {
    to: '/tools',
    label: 'Tools',
    icon: HiOutlineAdjustments,
  },
  {
    to: '/ocsr',
    label: 'OCSR',
    icon: HiOutlineCamera,
  },
  {
    to: '/about',
    label: 'About',
    icon: HiUserGroup,
    animated: true, // Flag to identify this icon for animation
  },
];

const Navigation = ({ isMobile = false, closeMenu = () => { } }) => {
  // CSS classes for different states of the navigation links
  const getLinkClasses = ({ isActive }) => {
    // Base classes for layout, padding, font, and transitions
    const baseClasses = `flex items-center rounded-md text-sm font-medium transition-all duration-200 ease-in-out group ${
      isMobile ? 'py-2 px-3 text-base' : 'px-3 py-2' // Slightly larger text/padding for mobile clarity
    }`;

    // Classes for the ACTIVE link state (Light & Dark)
    const activeClasses = isActive
      ? 'bg-sky-100 dark:bg-slate-700 text-sky-700 dark:text-white shadow-inner dark:shadow-none' // Example using theme colors
      : // Classes for the INACTIVE link state (Light & Dark)
      'text-slate-700 dark:text-slate-300 hover:bg-slate-200/60 dark:hover:bg-slate-700/60 hover:text-slate-900 dark:hover:text-white';

    return `${baseClasses} ${activeClasses}`;
  };

  // Add the animation keyframes to the document
  React.useEffect(() => {
    // Create a style element for our custom animations
    const style = document.createElement('style');
    style.innerHTML = `
      @keyframes usergroup-strobe {
        0%, 100% { opacity: 1; }
        50% { opacity: 0.5; }
      }
      
      @keyframes usergroup-shake {
        0%, 92%, 100% { transform: translateX(0); }
        94% { transform: translateX(-1px) rotate(-0.5deg); }
        96% { transform: translateX(1px) rotate(0.5deg); }
        98% { transform: translateX(-1px); }
      }
      
      .animate-usergroup-strobe {
        animation: usergroup-strobe 3s ease-in-out infinite, usergroup-shake 17s ease-in-out infinite;
      }
    `;
    // Append the style to the head
    document.head.appendChild(style);
    
    // Clean up on unmount
    return () => {
      document.head.removeChild(style);
    };
  }, []);

  return (
    // Adjusted spacing for desktop
    <nav className={`${isMobile ? 'flex flex-col space-y-1' : 'flex items-center space-x-1 lg:space-x-2'}`}>
      {navLinks.map((link) => {
        const Icon = link.icon; // Get the icon component
        const iconClassName = `flex-shrink-0 h-5 w-5 ${isMobile ? 'mr-3' : 'mr-1.5'}`;
        
        return (
          <NavLink
            key={link.to}
            to={link.to}
            className={getLinkClasses}
            onClick={isMobile ? closeMenu : undefined}
            end={link.exact} // Use end prop from config
          >
            {/* Render the animated StrobingIcon component for UserGroup icon, regular Icon for others */}
            {link.animated ? (
              <StrobingIcon icon={Icon} className={iconClassName} />
            ) : (
              <Icon className={iconClassName} aria-hidden="true" />
            )}
            <span>{link.label}</span>
          </NavLink>
        );
      })}
    </nav>
  );
};

export default Navigation;