import React from 'react'
import PropTypes from 'prop-types'

/**
 * For embedding YouTube videos.
 */
const YouTubeEmbed = ({ embedId }) => (
  <div className="video-responsive">
    <iframe
      title= 'Youtube iframe'
      width="607"
      height="340"
      src={`https://www.youtube.com/embed/${embedId}`}
      frameBorder="0"
      allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
      allowFullScreen
    />
  </div>
)

YouTubeEmbed.propTypes = {
  embedId: PropTypes.string.isRequired
}

export default YouTubeEmbed
