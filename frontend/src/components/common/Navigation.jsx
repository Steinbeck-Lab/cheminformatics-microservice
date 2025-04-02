// Description: This component renders a navigation menu with links to different sections of the application.
import React from 'react';
import { NavLink } from 'react-router-dom';
import {
  HiOutlineHome,
  HiOutlineBeaker,
  HiOutlineRefresh, // User's chosen icon
  HiOutlinePresentationChartLine, // User's chosen icon
  HiOutlineAdjustments, // User's chosen icon
  HiOutlineCamera
} from 'react-icons/hi';

// Navigation links configuration (using user's icons)
const navLinks = [
  {
    to: '/',
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
    icon: HiOutlineRefresh, // User's chosen icon
  },
  {
    to: '/depict',
    label: 'Depiction',
    icon: HiOutlinePresentationChartLine, // User's chosen icon
  },
  {
    to: '/tools',
    label: 'Tools',
    icon: HiOutlineAdjustments, // User's chosen icon
  },
  {
    to: '/ocsr',
    label: 'OCSR',
    icon: HiOutlineCamera,
  },
];

const Navigation = ({ isMobile = false, closeMenu = () => { } }) => {
  // CSS classes for different states of the navigation links
  const getLinkClasses = ({ isActive }) => {
    // Base classes for layout, padding, font, and transitions
    const baseClasses = `flex items-center rounded-md text-sm font-medium transition-all duration-200 ease-in-out group ${isMobile ? 'py-2 px-3 text-base' : 'px-3 py-2' // Slightly larger text/padding for mobile clarity
      }`;

    // Classes for the ACTIVE link state (Light & Dark)
    const activeClasses = isActive
      ? 'bg-sky-100 dark:bg-slate-700 text-sky-700 dark:text-white shadow-inner dark:shadow-none' // Example using theme colors
      : // Classes for the INACTIVE link state (Light & Dark)
      'text-slate-700 dark:text-slate-300 hover:bg-slate-200/60 dark:hover:bg-slate-700/60 hover:text-slate-900 dark:hover:text-white';

    return `${baseClasses} ${activeClasses}`;
  };

  return (
    // Adjusted spacing for desktop
    <nav className={`${isMobile ? 'flex flex-col space-y-1' : 'flex items-center space-x-1 lg:space-x-2'}`}>
      {navLinks.map((link) => {
        const Icon = link.icon; // Get the icon component
        return (
          <NavLink
            key={link.to}
            to={link.to}
            className={getLinkClasses}
            onClick={isMobile ? closeMenu : undefined}
            end={link.exact} // Use end prop from config
          >
            {/* Icon: Color should inherit from parent text color */}
            <Icon className={`flex-shrink-0 h-5 w-5 ${isMobile ? 'mr-3' : 'mr-1.5'}`} aria-hidden="true" />
            <span>{link.label}</span>
          </NavLink>
        );
      })}
    </nav>
  );
};

export default Navigation;
